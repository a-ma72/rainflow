#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <vector>
//#define RFC_MEM_ALLOC NULL
#define RFC_TP_STORAGE std::vector<RF::rfc_value_tuple_s>
#include <rainflow.hpp>



// Convert RFC error numbers to strings
static
const char* rfc_err_str( int nr )
{
    switch( nr )
    {
        case Rainflow::RFC_ERROR_NOERROR:            return "No error";
        case Rainflow::RFC_ERROR_INVARG:             return "Invalid arguments passed";
        case Rainflow::RFC_ERROR_UNSUPPORTED:        return "Unsupported feature";
        case Rainflow::RFC_ERROR_MEMORY:             return "Error on memory allocation";
        case Rainflow::RFC_ERROR_TP:                 return "Error while processing turning points";
        case Rainflow::RFC_ERROR_AT:                 return "Error while amplitude transformation";
        case Rainflow::RFC_ERROR_DH_BAD_STREAM:      return "Input stream must be unique";
        case Rainflow::RFC_ERROR_DH:                 return "Error while damage history calculation/access";
        case Rainflow::RFC_ERROR_LUT:                return "Error while accessing look up tables";
        case Rainflow::RFC_ERROR_DATA_OUT_OF_RANGE:  return "Input data leaves classrange";
        default:                                     return "Unexpected error";
    }
}


// Parse RFC counting parameters
static
int parse_rfc_kwargs( PyObject* kwargs, Py_ssize_t len, Rainflow *rf, Rainflow::rfc_res_method *res_method )
{
    PyObject   *empty           =  PyTuple_New(0);
    int         class_count     =  100;
    double      class_width     = -1;  // -1 = "calculated"
    double      class_offset    =  0;
    double      hysteresis      = -1;  // -1 = "calculated"
    int         enforce_margin  =  1;  // true
    int         use_hcm         =  0;  // false
    int         use_astm        =  0;  // false
    int         lc_method       =  0;  // Count rising slopes only
    int         flags           =  Rainflow::RFC_FLAGS_DEFAULT;
    int         auto_resize     =  0;  // false
    int         spread_damage   =  Rainflow::RFC_SD_TRANSIENT_23c;
    PyObject   *wl              =  NULL;
    double      wl_sd           =  1e3, wl_nd = 1e7, 
                wl_k            =  5,   wl_k2 = 5;

    *res_method = Rainflow::RFC_RES_REPEATED;

    char* kw[] = {"class_width", "class_count", "class_offset", 
                  "hysteresis","residual_method", "enforce_margin", "auto_resize",
                  "use_HCM", "use_ASTM", "spread_damage", "lc_method", "wl", NULL};

    if( !PyArg_ParseTupleAndKeywords( empty, kwargs, "d|iddiiiiiiiO", kw,
                                      &class_width,
                                      &class_count,
                                      &class_offset,
                                      &hysteresis,
                                       res_method,
                                      &enforce_margin,
                                      &auto_resize,
                                      &use_hcm,
                                      &use_astm,
                                      &spread_damage,
                                      &lc_method,
                                      &wl ) )
    {
        Py_DECREF( empty );
        return 0;
    }
    Py_DECREF( empty );

    // Parameters of the SN curve, if defined
    if( wl )
    {
        if( !PyDict_Check( wl ) )
        {
            PyErr_SetString( PyExc_RuntimeError, "Parameter 'wl' must be of type dict!" );
            return 0;
        }

        PyObject *key, *value;
        Py_ssize_t pos = 0;

        // Iterate over keys
        while( PyDict_Next( wl, &pos, &key, &value ) )
        {
            if( !PyUnicode_Check( key ) )
            {
                PyErr_SetString( PyExc_RuntimeError, "Only string keys allowed in wl dict!" );
                return 0;
            }

            if( PyUnicode_CompareWithASCIIString( key, "sd") == 0 )
            {
                wl_sd = PyFloat_AsDouble( value );
            }
            else if( PyUnicode_CompareWithASCIIString( key, "nd" ) == 0 )
            {
                wl_nd = PyFloat_AsDouble( value );
            }
            else if( PyUnicode_CompareWithASCIIString( key, "k" ) == 0 )
            {
                wl_k = PyFloat_AsDouble( value );
            }
            else if( PyUnicode_CompareWithASCIIString( key, "k2" ) == 0 )
            {
                wl_k2 = fabs( PyFloat_AsDouble( value ) );
            }
            else
            {
                PyErr_Format( PyExc_RuntimeError, "Wrong key used in wl dict: %O", key );
                return 0;
            }
        }

        if( wl_k2 < 0 ) wl_k2 = wl_k;
        if( hysteresis < 0 ) hysteresis = class_width;
    }

    if( !rf->init( class_count, class_width, class_offset, hysteresis, (Rainflow::rfc_flags_e)flags ) )
    {
        PyErr_Format( PyExc_RuntimeError, "Rainflow initialization error (%s)", rfc_err_str( rf->error_get() ) );
        return 0;
    }
    else
    {
        rf->flags_get( &flags );

        // lc_method
        flags &= ~Rainflow::RFC_FLAGS_COUNT_LC;

        switch( lc_method )
        {
            case 0:
                flags |= Rainflow::RFC_FLAGS_COUNT_LC_UP;
                break;

            case 1:
                flags |= Rainflow::RFC_FLAGS_COUNT_LC_DN;
                break;

            case 2:
                flags |= Rainflow::RFC_FLAGS_COUNT_LC;
                break;

            default:
                PyErr_SetString( PyExc_RuntimeError, "Parameter 'lc_method' must be 0, 1 or 2!" );
                return 0;
        }

        switch( auto_resize )
        {
            case 0:
                break;

            case 1:
                flags &= ~Rainflow::RFC_FLAGS_AUTORESIZE;
                flags |=  Rainflow::RFC_FLAGS_AUTORESIZE;
                break;

            default:
                PyErr_SetString( PyExc_RuntimeError, "Parameter 'auto_resize' must be 0 or 1!" );
                return 0;
        }

        rf->flags_set( flags, /* debugging */ false, /* overwrite */ true );
    }

    if( !rf->wl_init_modified( wl_sd, wl_nd, wl_k, wl_k2 ) )
    {
        PyErr_Format( PyExc_RuntimeError, "Rainflow initialization error (%s)", rfc_err_str( rf->error_get() ) );
        return 0;
    }

    if( spread_damage < (int)Rainflow::RFC_SD_NONE || spread_damage >= (int)Rainflow::RFC_SD_COUNT )
    {
        PyErr_SetString( PyExc_RuntimeError, "Unknown method for handling damage history!" );
        return 0;
    }

    if( spread_damage > (int)Rainflow::RFC_SD_NONE )
    {
        if( !rf->dh_init( (Rainflow::rfc_sd_method_e) spread_damage, NULL, (size_t)len, /*is_static*/ false ) )
        {
            PyErr_SetString( PyExc_MemoryError, "Error allocation damage history!" );
            return 0;
        }
    }

    if( (int)*res_method < (int)Rainflow::RFC_RES_NONE || (int)*res_method >= (int)Rainflow::RFC_RES_COUNT )
    {
        PyErr_SetString( PyExc_RuntimeError, "Unknown method for handling residue!" );
        return 0;
    }

    if( use_hcm && use_astm )
    {
        return 0;
    }

    if( use_hcm )
    {
        rf->ctx_get().counting_method = RF::RFC_COUNTING_METHOD_HCM;
    }

    if( use_astm )
    {
        rf->ctx_get().counting_method = RF::RFC_COUNTING_METHOD_ASTM;
    }

    return 1;
}


// Parse input data array
static 
int parse_rfc_input_series( PyObject* input_series_arg, PyArrayObject **arr_data, npy_double **data, Py_ssize_t *len )
{
    PyObject *arg1;

    *arr_data = NULL;
    *data = NULL;
    *len = 0;

    if( !PyArg_ParseTuple( input_series_arg, "O", &arg1 ) )
    {
        return 0;
    }
    /*
    if( !PyArray_Check( arg1 ) )
    {
        PyErr_SetString( PyExc_RuntimeError, "Not an ndarray!" );
        return 0;
    }
    */
    *arr_data = (PyArrayObject*)PyArray_FROM_OTF( arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
    if( *arr_data == NULL )
    {
        return 0;
    }

    // Number of dimensions
    int arr_nd = PyArray_NDIM( *arr_data );
    int arr_type = PyArray_TYPE( *arr_data );
    npy_intp *arr_dims = PyArray_DIMS( *arr_data );

    if( arr_nd != 1 )
    {
        PyErr_SetString( PyExc_RuntimeError, "Data must have only one dimension!" );
        return 0;
    }

    int r = PyArray_AsCArray( (PyObject**)arr_data, (void**)data, arr_dims, arr_nd, PyArray_DescrFromType(arr_type) );
    if( r < 0 )
    {
        PyErr_SetString( PyExc_RuntimeError, "Could not convert input to C array" );
        return 0;
    }

    *len = arr_dims[0];

    return 1;
}


// Process rainflow counting
static
int do_rainflow( Rainflow *rf, npy_double *data, Py_ssize_t len, Rainflow::rfc_res_method res_method )
{
    if( !rf->feed( data, len ) || !rf->finalize( res_method ) ) goto fail;

    return 1;
fail:
    PyErr_Format( PyExc_RuntimeError, "Error while counting (%s)", rfc_err_str( rf->error_get() ) );
    return 0;
}


// Prepare results
static
int prepare_results( Rainflow *rf, Rainflow::rfc_res_method res_method, PyObject **ret )
{
    const Rainflow::rfc_value_tuple_s *p_residue;
    Rainflow::rfc_counts_v ct;
    Rainflow::rfc_value_v sa;
    Rainflow::rfc_tp_storage tp;
    Rainflow::rfc_rfm_item_v rfm;
    unsigned u, class_count;
    double damage;
    const double *dh;
    size_t dh_cnt;
    PyArrayObject *arr;
    npy_intp len[2];

    *ret = NULL;

    // Retrieve range pair counts
    if( !rf->class_count( &u ) ) goto fail;
    class_count = u;
    if( !rf->rp_get( ct, sa ) )
    {
        goto fail_rfc;
    }

    // Create dict (return value)
    *ret = PyDict_New();
    if( *ret == NULL ) goto fail_cont;

    // Insert damage value
    if( !rf->damage( &damage ) ) goto fail_rfc;
    PyDict_SetItemString( *ret, "damage", PyFloat_FromDouble( damage ) );

    // Insert range counts
    len[0] = class_count;
    len[1] = 2;
    arr = (PyArrayObject*)PyArray_SimpleNew( 2, len, NPY_DOUBLE );
    if( !arr ) goto fail_cont;
    PyArray_FILLWBYTE( arr, 0 );
    for( unsigned i = 0; i < class_count; i++ )
    {
        *(double*)PyArray_GETPTR2( arr, i, 0 ) = (double)sa[i] * 2;  // range = 2 * amplitude
        *(double*)PyArray_GETPTR2( arr, i, 1 ) = (double)ct[i];
    }
    PyDict_SetItemString( *ret, "rp", (PyObject*)arr );
    Py_DECREF( arr );

    // Insert level crossings
    if( !rf->lc_get( ct, sa ) ) goto fail_rfc;
    len[0] = class_count;
    len[1] = 2;
    arr = (PyArrayObject*)PyArray_SimpleNew( 2, len, NPY_DOUBLE );
    if( !arr ) goto fail_cont;
    PyArray_FILLWBYTE( arr, 0 );
    for( unsigned i = 0; i < class_count; i++ )
    {
        *(double*)PyArray_GETPTR2( arr, i, 0 ) = (double)sa[i];  // class upper limit
        *(double*)PyArray_GETPTR2( arr, i, 1 ) = (double)ct[i];
    }
    PyDict_SetItemString( *ret, "lc", (PyObject*)arr );
    Py_DECREF( arr );

    // Insert turning points
    len[0] = rf->tp_storage().size();
    len[1] = 3;
    arr = (PyArrayObject*)PyArray_SimpleNew( 2, len, NPY_DOUBLE );
    if( !arr ) goto fail_cont;
    PyArray_FILLWBYTE( arr, 0 );
    for( size_t i = 0; i < rf->tp_storage().size(); i++ )
    {
        *(double*)PyArray_GETPTR2( arr, i, 0 ) = (double)rf->tp_storage()[i].pos;
        *(double*)PyArray_GETPTR2( arr, i, 1 ) = (double)rf->tp_storage()[i].value;
        *(double*)PyArray_GETPTR2( arr, i, 2 ) = (double)rf->tp_storage()[i].damage;
    }
    PyDict_SetItemString( *ret, "tp", (PyObject*)arr );
    Py_DECREF( arr );

    // Insert residue
    if( !rf->res_get( &p_residue, &u ) ) goto fail;
    len[0] = u;
    len[1] = 0;
    arr = (PyArrayObject*)PyArray_SimpleNew( 1, len, NPY_DOUBLE );
    if( !arr ) goto fail_cont;
    PyArray_FILLWBYTE( arr, 0 );
    for( unsigned i = 0; i < u; i++ )
    {
        *(double*)PyArray_GETPTR1( arr, i ) = (double)p_residue[i].value;
    }
    PyDict_SetItemString( *ret, "res", (PyObject*)arr );
    Py_DECREF( arr );

    // Insert rainflow matrix
    if( !rf->rfm_get( rfm ) ) goto fail_rfc;
    len[0] = class_count;
    len[1] = class_count;
    arr = (PyArrayObject*)PyArray_SimpleNew( 2, len, NPY_DOUBLE );
    if( !arr ) goto fail_cont;
    PyArray_FILLWBYTE( arr, 0 );
    for( size_t k = 0; k < rfm.size(); k++ )
    {
        unsigned i, j;
        i = rfm[k].from;
        j = rfm[k].to;
        *(double*)PyArray_GETPTR2( arr, i, j ) += (double)rfm[k].counts / RFC_FULL_CYCLE_INCREMENT;
    }
    PyDict_SetItemString( *ret, "rfm", (PyObject*)arr );
    Py_DECREF( arr );

    // Insert damage history
    if( !rf->dh_get( &dh, &dh_cnt ) ) goto fail_rfc;
    len[0] = dh_cnt;
    len[1] = 0;
    arr = (PyArrayObject*)PyArray_SimpleNew( 1, len, NPY_DOUBLE );
    if( !arr ) goto fail_cont;
    PyArray_FILLWBYTE( arr, 0 );
    for( size_t i = 0; i < dh_cnt; i++ )
    {
        *(double*)PyArray_GETPTR1( arr, i ) = dh[i];
    }
    PyDict_SetItemString( *ret, "dh", (PyObject*)arr );
    Py_DECREF( arr );

    return 1;

fail:
    PyErr_SetString( PyExc_MemoryError, "Preparing range pair counting" );
    goto fail_cont;
fail_rfc:
    PyErr_Format( PyExc_MemoryError, "Preparing range pair counting (%s)", rfc_err_str( rf->error_get() ) );
    goto fail_cont;
fail_cont:
    if( *ret ) Py_DECREF( *ret );
    return 0;
}


static PyObject* rfc( PyObject *self, PyObject *args, PyObject *kwargs )
{
    PyArrayObject *arr_data = NULL;
    PyObject *ret = NULL;
    npy_double *data = NULL;
    Rainflow rf;
    Rainflow::rfc_res_method res_method;
    Py_ssize_t len;
    bool ok = false;

    do
    {
        if( !parse_rfc_input_series( args, &arr_data, &data, &len ) )
        {
            break;
        }

        if( !parse_rfc_kwargs( kwargs, len, &rf, &res_method ) )
        {
            break;
        }

        if( !do_rainflow( &rf, data, len, res_method ) )
        {
            break;
        }

        if( !prepare_results( &rf, res_method, &ret ) )
        {
            break;
        }
        
        ok = true;
    }
    while(0);


    if( !ok && ret )
    {
        Py_DECREF( ret );
        ret = NULL;
    }
    
    rf.deinit();

    if( arr_data && data )
    {
        Py_DECREF( arr_data );
        PyArray_Free( (PyObject*)arr_data, (void*)data );
    }

    return ret;
}


// Exported methods are collected in a table
PyMethodDef method_table[] = {
    {"rfc", (PyCFunction) rfc, METH_VARARGS | METH_KEYWORDS, "Rainflow counting"},
    {NULL, NULL, 0, NULL} // Sentinel value ending the table
};


// A struct contains the definition of a module
PyModuleDef mymath_module = {
    PyModuleDef_HEAD_INIT,
    "rfcnt", // Module name
    "Rainflow counting module",
    -1,   // Optional size of the module state memory
    method_table,
    NULL, // Optional slot definitions
    NULL, // Optional traversal function
    NULL, // Optional clear function
    NULL  // Optional module deallocation function
};


// The module init function
PyMODINIT_FUNC PyInit_rfcnt(void) {
    PyObject* mod = PyModule_Create(&mymath_module);
    // Initialize numpy
    import_array();
    return mod;
}
