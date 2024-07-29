#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
// C ABI compatibility: https://stackoverflow.com/a/74296491/11492317
// https://pypi.org/project/oldest-supported-numpy/
#include <numpy/arrayobject.h>
// #include "numpy_arrayobject.h"
#include <vector>
//#define RFC_MEM_ALLOC nullptr
#define RFC_TP_STORAGE std::vector<RF::rfc_value_tuple_s>
#include <rainflow.hpp>

typedef std::vector<Rainflow::rfc_value_tuple_s> rfc_residuum_vec;




static
PyObject* _numpy_api_version(PyObject *self )
{
    // Return the NumPy API version as a Python integer object
    return Py_BuildValue("i", (int)(NPY_API_VERSION));
}


static
bool convert_to_numpy_array(PyObject* input_series_arg, PyArrayObject **arr_data, npy_double **data, Py_ssize_t *len )
{
    *arr_data = nullptr;
    *data = nullptr;
    *len = 0;

    /*
    if( !PyArray_Check( arg1 ) )
    {
        PyErr_SetString( PyExc_RuntimeError, "Not an ndarray!" );
        return 0;
    }
    */
    *arr_data = (PyArrayObject*)PyArray_FROM_OTF( input_series_arg, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
    if( *arr_data == nullptr )
    {
        return false;
    }

    // Number of dimensions
    int arr_nd = PyArray_NDIM( *arr_data );
    int arr_type = PyArray_TYPE( *arr_data );
    npy_intp *arr_dims = PyArray_DIMS( *arr_data );

    if( arr_nd != 1 )
    {
        PyErr_SetString( PyExc_RuntimeError, "Data must have only one dimension!" );
        return false;
    }

    int r = PyArray_AsCArray( (PyObject**)arr_data, (void**)data, arr_dims, arr_nd, PyArray_DescrFromType(arr_type) );
    if( r < 0 )
    {
        PyErr_SetString( PyExc_RuntimeError, "Could not convert input to C array" );
        return false;
    }

    *len = arr_dims[0];

    return true;
}


static
bool get_numeric_double(PyObject *Py_value, const char *name, double &value )
{
    // Check if the value is a number
    if( PyFloat_Check(Py_value) || PyLong_Check(Py_value) )
    {
        value = PyFloat_AsDouble(Py_value);
        if( PyErr_Occurred() )
        {
            PyErr_Format( PyExc_ValueError, "Invalid value for `%s`", name );
            return false;
        }
    } else {
        PyErr_Format( PyExc_TypeError, "`%s` must be a numeric type", name );
        return false;
    }

    return true;
}


static
bool get_dict_item_double( PyObject *dict, const char* name, double &value, double default_value )
{
    // Check if the argument is a dictionary
    if( !PyDict_Check(dict) )
    {
        PyErr_Format( PyExc_TypeError, "Argument `%s` must be a dictionary.", name );
        return false;
    }

    PyObject *Py_value = PyDict_GetItemString( dict, name );
    if( Py_value != nullptr )
    {
        if( !get_numeric_double(Py_value, name, value) ) return false;
    }
    else
    {
        value = default_value;
    }

    return true;
}


static
bool get_dict_wl( PyObject *Py_wl, const char *name, Rainflow::rfc_wl_param_s &wl )
{
    if( !PyDict_Check( Py_wl ) )
    {
        PyErr_Format( PyExc_RuntimeError, "Parameter '%' must be of type dict!", name );
        return false;
    }

    PyObject *key, *value;
    Py_ssize_t pos = 0;
    bool wl_k2_set = false;

    // Iterate over keys
    while( PyDict_Next( Py_wl, &pos, &key, &value ) )
    {
        if( !PyUnicode_Check( key ) )
        {
            PyErr_Format( PyExc_RuntimeError, "Only string keys allowed in `%s`", name );
            return false;
        }

        if( PyUnicode_CompareWithASCIIString( key, "sd") == 0 )
        {
            if( !get_numeric_double(value, "sd", wl.sd) ) return false;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "nd" ) == 0 )
        {
            if( !get_numeric_double(value, "nd", wl.nd) ) return false;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "k" ) == 0 )
        {
            if( !get_numeric_double(value, "k", wl.k) ) return false;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "k2" ) == 0 )
        {
            if( !get_numeric_double(value, "k2", wl.k2) ) return false;
            wl_k2_set = true;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "sx" ) == 0 )
        {
            if( !get_numeric_double(value, "sx", wl.sx) ) return false;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "nx" ) == 0 )
        {
            if( !get_numeric_double(value, "nx", wl.nx) ) return false;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "omission" ) == 0 )
        {
            if( !get_numeric_double(value, "omission", wl.omission) ) return false;
        }
        else
        {
            PyErr_Format( PyExc_RuntimeError, "Wrong key used in wl dict: `%S`", key );
            return false;
        }
    }

    if( !wl_k2_set )
    {
        wl.k2 = wl.k;
    }
    return true;
}


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
        case Rainflow::RFC_ERROR_DATA_OUT_OF_RANGE:  return "Input data leaves class range";
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
    PyObject   *wl              =  nullptr;
    double      wl_sd           =  1e3, wl_nd = 1e7,
                wl_k            =  5,   wl_k2 = 5;

    *res_method = Rainflow::RFC_RES_REPEATED;

    const char* kw[] = {"class_width", "class_count", "class_offset",
                        "hysteresis","residual_method", "enforce_margin", "auto_resize",
                        "use_HCM", "use_ASTM", "spread_damage", "lc_method", "wl", nullptr};

    if( !PyArg_ParseTupleAndKeywords( empty, kwargs, "d|iddi$ppppiiO", (char**)kw,
                                      &class_width,     // d
                                      &class_count,     // i
                                      &class_offset,    // d
                                      &hysteresis,      // d
                                       res_method,      // i
                                      &enforce_margin,  // p
                                      &auto_resize,     // p
                                      &use_hcm,         // p
                                      &use_astm,        // p
                                      &spread_damage,   // i
                                      &lc_method,       // i
                                      &wl ) )           // O
    {
        Py_DECREF( empty );
        return 0;
    }
    Py_DECREF( empty );

    // Parameters of the SN curve, if defined
    if( wl && wl != Py_None )
    {
        Rainflow::rfc_wl_param_s wl_param = {0};
        if( get_dict_wl( wl, "wl", wl_param) )
        {
            wl_sd = wl_param.sd;
            wl_nd = wl_param.nd;
            wl_k = wl_param.k;
            wl_k2 = wl_param.k2;

            if( wl_param.sx != 0.0 || wl_param.nx != 0.0 || wl_param.omission != 0.0 )
            {
                PyErr_Warn(PyExc_Warning, "Keys `sx`, `nx` and `omission` get ignored.");
            }
        }
        else
        {
            return 0;
        }
    }

    if( hysteresis < 0 ) hysteresis = class_width;

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
                flags &= ~Rainflow::RFC_FLAGS_AUTORESIZE;
                break;

            case 1:
                flags |=  Rainflow::RFC_FLAGS_AUTORESIZE;
                break;

            default:
                PyErr_SetString( PyExc_RuntimeError, "Parameter 'auto_resize' must be 0 or 1!" );
                return 0;
        }

        switch( enforce_margin )
        {
            case 0:
                flags &= ~Rainflow::RFC_FLAGS_ENFORCE_MARGIN;
                break;

            case 1:
                flags |=  Rainflow::RFC_FLAGS_ENFORCE_MARGIN;
                break;

            default:
                PyErr_SetString( PyExc_RuntimeError, "Parameter 'enforce_margin' must be 0 or 1!" );
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
        if( !rf->dh_init( (Rainflow::rfc_sd_method_e) spread_damage, nullptr, (size_t)len, /*is_static*/ false ) )
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


// Process rainflow counting
static
int do_rainflow( Rainflow *rf, npy_double *data, Py_ssize_t len, Rainflow::rfc_res_method res_method, rfc_residuum_vec &residuum_raw )
{
    const Rainflow::rfc_value_tuple_s *residuum;
    unsigned residuum_len;

    if( !rf->feed( data, len ) ) goto fail;

    if( !rf->res_get( &residuum, &residuum_len ) ) goto fail;

    residuum_raw = rfc_residuum_vec( residuum, residuum + residuum_len );

    // (With regard to finalize_res_repeated() in rainflow.c, remove pending cycle)
    if( residuum_raw.size() >= 4 )
    {
        size_t idx = residuum_raw.size() - 4;

        unsigned A = residuum_raw[idx+0].cls;
        unsigned B = residuum_raw[idx+1].cls;
        unsigned C = residuum_raw[idx+2].cls;
        unsigned D = residuum_raw[idx+3].cls;

        if( B > C )
        {
            unsigned temp = B;
            B = C;
            C = temp;
        }

        if( A > D )
        {
            unsigned temp = A;
            A = D;
            D = temp;
        }

        /* Check for closed cycles [3] */
        if( A <= B && C <= D )
        {
            // Remove points B and C
            residuum_raw.erase(residuum_raw.begin() + 1, residuum_raw.begin() + 3);
        }
    }

    if( !rf->finalize( res_method ) ) goto fail;

    return 1;
fail:
    PyErr_Format( PyExc_RuntimeError, "Error while counting (%s)", rfc_err_str( rf->error_get() ) );
    return 0;
}


// Prepare results
static
int prepare_results( Rainflow *rf, Py_ssize_t data_len, Rainflow::rfc_res_method res_method, rfc_residuum_vec &residuum_raw, PyObject **ret )
{
    const Rainflow::rfc_value_tuple_s *p_residue;
    Rainflow::rfc_counts_v ct;
    Rainflow::rfc_value_v sa;
    Rainflow::rfc_tp_storage tp;
    Rainflow::rfc_rfm_item_v rfm;
    Rainflow::rfc_wl_param_s wl;
    unsigned u, class_count;
    double damage;
    const double *dh;
    size_t dh_cnt;
    PyArrayObject *arr;
    PyObject *obj;
    npy_intp len[2];

    *ret = nullptr;

    // Retrieve range pair counts
    if( !rf->class_count( &u ) ) goto fail;
    class_count = u;
    if( !rf->rp_get( ct, sa ) )
    {
        goto fail_rfc;
    }

    // Create dict (return value)
    *ret = PyDict_New();
    if( *ret == nullptr ) goto fail_cont;

    // Insert damage value
    if( !rf->damage( &damage ) ) goto fail_rfc;
    PyDict_SetItemString( *ret, "damage", PyFloat_FromDouble( damage ) );

    // Insert range pair counts
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

    // Insert residue_raw
    u = (unsigned int)residuum_raw.size();
    len[0] = u;
    len[1] = 0;
    arr = (PyArrayObject*)PyArray_SimpleNew( 1, len, NPY_DOUBLE );
    if( !arr ) goto fail_cont;
    PyArray_FILLWBYTE( arr, 0 );
    for( unsigned i = 0; i < u; i++ )
    {
        *(double*)PyArray_GETPTR1( arr, i ) = (double)residuum_raw[i].value;
    }
    PyDict_SetItemString( *ret, "res_raw", (PyObject*)arr );
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
    len[0] = data_len;
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

    // WL impaired
    if( !rf->wl_param_get_impaired( wl ) ) goto fail_rfc;
    obj = PyDict_New();
    if( !obj ) goto fail_rfc;
    PyDict_SetItemString( obj, "sx", PyFloat_FromDouble( wl.sx ) );
    PyDict_SetItemString( obj, "nx", PyFloat_FromDouble( wl.nx ) );
    PyDict_SetItemString( obj, "sd", PyFloat_FromDouble( wl.sd ) );
    PyDict_SetItemString( obj, "nd", PyFloat_FromDouble( wl.nd ) );
    PyDict_SetItemString( obj, "k", PyFloat_FromDouble( wl.k ) );
    PyDict_SetItemString( obj, "k2", PyFloat_FromDouble( wl.k2 ) );
    PyDict_SetItemString( obj, "q", PyFloat_FromDouble( wl.q ) );
    PyDict_SetItemString( obj, "q2", PyFloat_FromDouble( wl.q2 ) );
    PyDict_SetItemString( obj, "omission", PyFloat_FromDouble( wl.omission ) );
    PyDict_SetItemString( obj, "D", PyFloat_FromDouble( wl.D ) );
    PyDict_SetItemString( *ret, "wl_miner_consequent", obj );

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


static
PyObject* rfc( PyObject *self, PyObject *args, PyObject *kwargs )
{
    PyArrayObject *arr_data = nullptr;
    PyObject *arg1, *ret = nullptr;
    npy_double *data = nullptr;
    Rainflow rf;
    Rainflow::rfc_res_method res_method;
    rfc_residuum_vec residuum_raw;
    Py_ssize_t len;
    bool ok = false;

    if( !PyArg_ParseTuple( args, "O", &arg1 ) )
    {
        return nullptr;
    }

    do
    {
        if( !convert_to_numpy_array(arg1, &arr_data, &data, &len) )
        {
            break;
        }

        if( !parse_rfc_kwargs( kwargs, len, &rf, &res_method ) )
        {
            break;
        }

        if( !do_rainflow( &rf, data, len, res_method, residuum_raw ) )
        {
            break;
        }

        if( !prepare_results( &rf, len, res_method, residuum_raw, &ret ) )
        {
            break;
        }

        ok = true;
    }
    while(0);


    if( !ok && ret )
    {
        Py_DECREF( ret );
        ret = nullptr;
    }

    rf.deinit();

    if( arr_data && data )
    {
        Py_DECREF( arr_data );
        PyArray_Free( (PyObject*)arr_data, (void*)data );
    }

    return ret;
}


static
PyObject* damage_from_rp( PyObject *self, PyObject *args, PyObject *kwargs )
{
    PyObject                    *Py_Sa, *Py_counts,
                                *Py_wl = nullptr;
    PyArrayObject               *arr = nullptr;
    double                      *buffer;
    Py_ssize_t                  size;
    unsigned int                class_count = 0;
    Rainflow::rfc_double_v      vec_sa;
    Rainflow::rfc_counts_v      vec_counts;
    std::vector<npy_intp>       vec_sorted_indices;
    Rainflow::rfc_wl_param_s    wl = {0};
    enum Rainflow::rfc_rp_damage_method damage_method = Rainflow::RFC_RP_DAMAGE_CALC_METHOD_DEFAULT;
    const char* kw[] = { "Sa", "counts",
                         "wl",
                         "method",
                         nullptr};

    if( !PyArg_ParseTupleAndKeywords( args, kwargs, "OO|$Oi:damage_from_rp", (char**)kw,
                                      &Py_Sa, &Py_counts, &Py_wl, &damage_method ) )
    {
        return nullptr;
    }

    if( damage_method < 0 || damage_method > 3 )
    {
        PyErr_SetString(PyExc_ValueError, "´damage_method` must be in range 0 to 3.");
        return nullptr;
    }

    wl.sx = 1e3;
    wl.nx = 1e7;
    wl.k = wl.k2 = 5;
    if( Py_wl )
    {
        if( !get_dict_wl( Py_wl, "wl", wl ) ) return nullptr;
    }

    if(convert_to_numpy_array(Py_Sa, &arr, &buffer, &size) )
    {
        PyObject *sorted_indices = PyArray_ArgSort( arr, 0, NPY_QUICKSORT );
        if( !sorted_indices || (PyObject*)sorted_indices == Py_None )
        {
            PyErr_SetString( PyExc_RuntimeError, "Error while sorting `Sa`." );
        }
        else
        {
            npy_intp* sorted_indices_data = (npy_intp*) PyArray_DATA( (PyArrayObject*)sorted_indices );
            vec_sorted_indices.assign( sorted_indices_data, sorted_indices_data + size );
            for( auto i : vec_sorted_indices )
            {
                vec_sa.push_back( buffer[i] );
            }
        }
        Py_XDECREF( arr );
        Py_XDECREF( sorted_indices );
    }
    else
    {
        PyErr_SetString( PyExc_ValueError, "Unable to convert input array `Sa`.");
        Py_XDECREF( arr );
        return nullptr;
    }

    if(convert_to_numpy_array(Py_counts, &arr, &buffer, &size) )
    {
        if( vec_sa.size() != (size_t)size )
        {
            PyErr_SetString( PyExc_ValueError, "`Sa` and `counts` must be of same size." );
        }
        else {
            double value;
            class_count = size;
            for (auto i: vec_sorted_indices)
            {
                value = buffer[i];
                if (value < 0)
                {
                    PyErr_SetString(PyExc_ValueError, "Negative values in `counts`.");
                    break;
                }
                vec_counts.push_back( (Rainflow::rfc_counts_t)(value + 0.5) );
            }
            Py_XDECREF(arr);
        }
    }
    else
    {
        PyErr_SetString( PyExc_ValueError, "Unable to convert input array `Sa`.");
        Py_XDECREF( arr );
    }

    if( PyErr_Occurred() )
    {
        return nullptr;
    }

    do {
        Rainflow rf;
        double damage;

        if( !rf.init( class_count, 1, 0, 0 ) ) break;
        if( !rf.wl_init_any( &wl ) ) break;
        if( !rf.damage_from_rp( damage, vec_counts, vec_sa, damage_method ) ) break;

        rf.deinit();
        return Py_BuildValue( "f", damage );
    } while(0);

    PyErr_SetString( PyExc_RuntimeError, "Error while calculation.");
    return nullptr;
}


// Exported methods are collected in a table
PyMethodDef method_table[] = {
    {"_numpy_api_version", (PyCFunction) _numpy_api_version, METH_NOARGS, "NumPy API version"},
    {"rfc", (PyCFunction) rfc, METH_VARARGS | METH_KEYWORDS, "Rainflow counting"},
    {"damage_from_rp", (PyCFunction) damage_from_rp, METH_VARARGS | METH_KEYWORDS, "Damage accumulation from range pair counting."},
    {nullptr, nullptr, 0,                                                       nullptr} // Sentinel value ending the table
};


// A struct contains the definition of a module
PyModuleDef mymath_module = {
    PyModuleDef_HEAD_INIT,
    "rfcnt", // Module name
    "Rainflow counting module",
    -1,   // Optional size of the module state memory
    method_table,
    nullptr, // Optional slot definitions
    nullptr, // Optional traversal function
    nullptr, // Optional clear function
    nullptr  // Optional module deallocation function
};


// The module init function
PyMODINIT_FUNC PyInit_rfcnt(void) {
    PyObject* mod = PyModule_Create(&mymath_module);
    // Initialize numpy
    import_array();
    return mod;
}
