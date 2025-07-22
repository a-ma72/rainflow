#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
// C ABI compatibility: https://stackoverflow.com/a/74296491/11492317
// https://pypi.org/project/oldest-supported-numpy/
#include <numpy/arrayobject.h>
// #include "numpy_arrayobject.h"
#include <vector>
#include <string>
//#define RFC_MEM_ALLOC nullptr
#define RFC_TP_STORAGE std::vector<RF::rfc_value_tuple_s>
#include <rainflow.hpp>

typedef std::vector<Rainflow::rfc_value_tuple_s> rfc_residuum_vec;




/**
 * @brief Returns the NumPy API version.
 *
 * This function retrieves the current version of the NumPy C API
 * and returns it as a Python integer object.
 *
 * @param self The self-reference object, usually a placeholder for bound methods.
 *             It is unused in this static method.
 * @return A Python integer object representing the NumPy C API version.
 */
static
PyObject* _numpy_api_version(PyObject *self )
{
    // Return the NumPy API version as a Python integer object
    return Py_BuildValue("i", (int)(NPY_API_VERSION));
}


/**
 * Converts a Python object to a NumPy array and extracts its data.
 *
 * This function attempts to interpret a Python object as a 1-dimensional
 * NumPy array of double precision floats (`npy_double`). If successful, it
 * returns a pointer to the NumPy array and the underlying data in C array
 * format. Additionally, the length of the array is extracted and returned.
 *
 * @param input_series_arg A Python object to be interpreted as a NumPy array.
 *                         Must be compatible with a 1-dimensional array of double precision floats.
 * @param arr_data A pointer to a `PyArrayObject*` that will store the resulting NumPy array if conversion
 *                 is successful, or `nullptr` if not.
 * @param data A pointer to an `npy_double*` that will point to the underlying C array data of the NumPy
 *             array upon successful conversion.
 * @param len A pointer to a `Py_ssize_t` that will store the length of the array (number of elements)
 *            upon successful conversion.
 *
 * @return `true` if the Python object is successfully converted to a 1-dimensional double precision
 *         NumPy array and the data is extracted; `false` otherwise.
 *
 * @note If the conversion fails at any step (e.g., data is not 1-dimensional, incompatible
 *       types, or memory issues), a Python exception is raised, and the pointers `arr_data`,
 *       `data`, and `len` will be set to `nullptr` or 0.
 *
 * @warning The function does not perform a full argument validation at the beginning, so it is
 *          caller's responsibility to ensure the passed object is at least a valid Python object.
 */
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


/**
 * @brief Extracts a numeric value (float or int) from a Python object and converts it to a double.
 *
 * This function checks whether the given Python object represents a numeric value (either
 * `int` or `float`). If it does, the value is converted to a C++ `double` and assigned to the provided
 * reference. If the object is not numeric, or the conversion fails, an appropriate Python exception
 * is raised and `false` is returned.
 *
 * @param Py_value The Python object expected to be a numeric value.
 * @param name A C-string used for error messages to indicate the name of the parameter being checked.
 * @param value A reference to a `double` that will receive the converted value if successful.
 *
 * @return `true` if the conversion was successful; `false` otherwise.
 *
 * @throws Raises a Python `TypeError` if `Py_value` is not an int or float.
 *         Raises a `ValueError` if the conversion to `double` fails.
 */
static
bool get_numeric_double( PyObject *Py_value, const char *name, double &value )
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


/**
 * @brief Validates and adjusts a named parameter value.
 *
 * This function checks and potentially modifies a numerical value based on its associated name.
 * It is typically used to enforce constraints or normalize parameters in a Python extension module.
 *
 * @param value Reference to the numeric value to validate or modify.
 * @param name Reference to the name of the parameter.
 *             - For "sd", "nd", "sx", "nx" or "omission": value must be non-negative.
 *               If negative, sets a Python ValueError and returns false.
 *             - For "k" or "k2": value is replaced with its absolute value.
 *             - Other names are accepted without validation or modification.
 *
 * @return true if the value is valid or was successfully adjusted; false if validation fails.
 *
 * @note This function uses the Python C API (PyErr_Format) to set errors.
 *       It is marked static and intended for internal use only.
 */
static
bool check_named_value( double &value, const std::string &name )
{
    if( name == "sd" || name == "nd" || name == "sx" || name == "nx" || name == "omission" )
    {
        if( value < 0.0 )
        {
            PyErr_Format( PyExc_ValueError, "Invalid value for `%s`", name.c_str() );
            return false;
        }
    }
    else if( name == "k" || name == "k2" )
    {
        value = fabs( value );
    }

    return true;
}


/**
 * @brief Extracts a double value from a Python dictionary, or assigns a default if key is missing.
 *
 * This function attempts to retrieve a key-value pair from a given Python dictionary (`dict`) where
 * the value is expected to be numeric (either an integer or float). If the key exists and the value
 * is valid, it is converted to a `double` and assigned to the provided reference `value`. If the key
 * does not exist, the `default_value` is assigned instead. If the value exists but is not numeric,
 * a Python exception is raised.
 *
 * @param dict A Python dictionary object to query.
 * @param name The name of the key to look for in the dictionary.
 * @param value A reference to a `double` that will hold the extracted or default value.
 * @param default_value The value to assign to `value` if the key is not present in the dictionary.
 *
 * @return `true` if the value was successfully extracted or the default assigned;
 *         `false` if an error occurred (e.g. dictionary is invalid or value is of incorrect type).
 *
 * @throws Raises a Python `TypeError` if `dict` is not a dictionary, or if the value for `name` is not numeric.
 *         Raises a `ValueError` if the numeric value cannot be converted to a `double`.
 */
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


/**
 * @brief Parses a Python dictionary representing Wöhler (S-N curve) parameters into a C++ struct.
 *
 * This function validates and extracts key-value pairs from a Python dictionary intended to define
 * Wöhler curve parameters. Recognized keys include `"sd"`, `"nd"`, `"k"`, `"k2"`, `"sx"`, `"nx"`, and `"omission"`,
 * each of which must have a numeric value. The extracted values are stored in the `wl` output struct.
 *
 * If `"k2"` is not provided, it is automatically set to the value of `"k"`.
 * Unknown keys result in a runtime error.
 *
 * If `"Py_wl"` is not provided, `"wl"` is left untouched.
 *
 * @param Py_wl The Python dictionary containing Wöhler parameters or Py_None.
 * @param name A name used in error messages to identify the context of the parameter.
 * @param wl A reference to a `Rainflow::rfc_wl_param_s` struct where parsed values will be stored.
 *
 * @return `true` if all required and optional parameters were successfully parsed or set;
 *         `false` if the input is invalid, a required value is incorrect, or if an unsupported key is encountered.
 *
 * @throws Raises a Python `RuntimeError`, `TypeError`, or `ValueError` depending on the specific issue:
 *         - Non-dictionary input
 *         - Non-string keys
 *         - Missing or invalid numeric values
 *         - Unsupported keys in the dictionary
 */
static
bool get_dict_wl( PyObject *Py_wl, const char *name, Rainflow::rfc_wl_param_s &wl, bool *extended_def = nullptr )
{
    /* if( Py_IsNone( Py_wl ) ) */
    if( Py_wl == Py_None )
    {
        return true;
    }

    if( !PyDict_Check( Py_wl ) )
    {
        PyErr_Format( PyExc_RuntimeError, "Parameter '%s' must be of type dict!", name );
        return false;
    }

    if( extended_def )
    {
        *extended_def = false;
    }

    PyObject *key, *value;
    Py_ssize_t pos = 0;

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
            if( !get_numeric_double( value, "sd", wl.sd) ) return false;
            if( !check_named_value( wl.sd, "sd" ) ) return false;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "nd" ) == 0 )
        {
            if( !get_numeric_double(value, "nd", wl.nd) ) return false;
            if( !check_named_value( wl.nd, "nd" ) ) return false;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "k" ) == 0 )
        {
            if( !get_numeric_double(value, "k", wl.k) ) return false;
            if( !check_named_value( wl.k, "k" ) ) return false;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "k2" ) == 0 )
        {
            if( !get_numeric_double(value, "k2", wl.k2) ) return false;
            if( !check_named_value( wl.k2, "k2" ) ) return false;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "sx" ) == 0 )
        {
            if( !get_numeric_double(value, "sx", wl.sx) ) return false;
            if( !check_named_value( wl.sx, "sx" ) ) return false;
            if( extended_def ) *extended_def = true;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "nx" ) == 0 )
        {
            if( !get_numeric_double(value, "nx", wl.nx) ) return false;
            if( !check_named_value( wl.nx, "nx" ) ) return false;
            if( extended_def ) *extended_def = true;
        }
        else if( PyUnicode_CompareWithASCIIString( key, "omission" ) == 0 )
        {
            if( !get_numeric_double(value, "omission", wl.omission) ) return false;
            if( !check_named_value( wl.omission, "omission" ) ) return false;
            if( extended_def ) *extended_def = true;
        }
        else
        {
            PyErr_Format( PyExc_RuntimeError, "Wrong key used in wl dict: `%S`", key );
            return false;
        }
    }

    wl.k = fabs( wl.k );
    wl.k2 = fabs( wl.k2 );

    if( wl.k == 0.0 )
    {
        wl.k = 5;
    }

    if( wl.sd == 0.0 && wl.sx == 0.0 )
    {
        wl.sx = 1e3;
    }

    if( wl.nd == 0.0 && wl.nx == 0.0 )
    {
        wl.nx = 1e7;
    }

    if( wl.k2 < wl.k )
    {
        wl.k2 = wl.k;
    }

    if( wl.sx == 0.0 && wl.sd > 0.0 )
    {
        wl.sx = wl.sd;
        wl.sd = 0.0;
    }

    if( wl.nx == 0.0 && wl.nd > 0.0 )
    {
        wl.nx = wl.nd;
        wl.nd = DBL_MAX;
    }

    return true;
}


/**
 * @brief Converts a Rainflow error code into a human-readable string.
 *
 * This function maps predefined error codes from the Rainflow module
 * to descriptive error messages, which are returned as C-style strings.
 * It is typically used for diagnostic and error reporting purposes.
 *
 * @param nr The error code as defined in the `Rainflow::RFC_ERROR_*` enum.
 *
 * @return A constant string describing the error corresponding to the input code.
 *         If the code is unrecognized, a default message "Unexpected error" is returned.
 */
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


/**
 * @brief Parses keyword arguments for Rainflow counting configuration.
 *
 * This function reads and interprets parameters from a Python dictionary (`kwargs`)
 * to configure a `Rainflow` instance for cycle counting. It supports configuration
 * of class parameters, counting options, and Wöhler (S-N curve) settings. Optional
 * and required parameters are validated and converted to appropriate types.
 *
 * Key parameters (with defaults) include:
 * - `class_width` (required): Width of each counting class (float).
 * - `class_count` (default=100): Number of classes (int).
 * - `class_offset` (default=0.0): Offset applied to class boundaries (float).
 * - `hysteresis` (default=-1): Threshold below which cycles are ignored. If < 0, uses class_width (float).
 * - `residual_method` (default=RFC_RES_REPEATED): How residuals are handled (int enum).
 * - `enforce_margin` (default=1): Whether to enforce class range margins (bool).
 * - `auto_resize` (default=0): Whether to auto-expand the class range (bool).
 * - `use_HCM` / `use_ASTM` (default=0): Selects counting method (bool, mutually exclusive).
 * - `spread_damage` (default=RFC_SD_TRANSIENT_23c): Damage history method (int enum).
 * - `lc_method` (default=0): Level crossing method (0: up, 1: down, 2: both).
 * - `wl`: Dictionary defining Wöhler parameters (`sd`, `nd`, `k`, `k2`, etc.).
 *
 * If a Wöhler dictionary is provided, additional validation is performed. Optional
 * keys like `sx`, `nx`, and `omission` are accepted but ignored with a warning.
 *
 * @param kwargs The Python dictionary containing the keyword arguments.
 * @param len The length of the input time series data (used for damage history sizing).
 * @param rf A pointer to a `Rainflow` instance to be configured.
 * @param res_method A pointer to store the selected residual method.
 *
 * @return `1` on success; `0` if any validation, memory allocation, or logic fails.
 *
 * @throws Raises various Python exceptions (`TypeError`, `RuntimeError`, etc.)
 *         with descriptive messages when validation or initialization fails.
 */
static
int parse_rfc_kwargs( PyObject* kwargs, Py_ssize_t len, Rainflow *rf, Rainflow::rfc_res_method *res_method )
{
    PyObject   *empty           =  PyTuple_New(0);
    int         class_count     =  100;
    double      class_width     =  1;  // will be overwritten by required argument
    double      class_offset    =  0;
    double      hysteresis      = -1;  // -1 => "use class_width as hysteresis"
    int         enforce_margin  =  1;  // true
    int         use_hcm         =  0;  // false
    int         use_astm        =  0;  // false
    int         lc_method       =  0;  // Count rising slopes only
    int         flags           =  Rainflow::RFC_FLAGS_DEFAULT;
    int         auto_resize     =  0;  // false
    int         spread_damage   =  Rainflow::RFC_SD_TRANSIENT_23c;
    PyObject   *wl              =  nullptr;
    double      wl_sx           =  1e3, wl_nx = 1e7,
                wl_sd           =  0.0, wl_nd = DBL_MAX,
                wl_k            =  5,   wl_k2 = 5,
                wl_omission     =  0.0;
    bool        wl_extended_def =  false;

    *res_method = Rainflow::RFC_RES_REPEATED;

    const char* kw[] = {"class_width", "class_count", "class_offset",
                        "hysteresis","residual_method", "enforce_margin", "auto_resize",
                        "use_HCM", "use_ASTM", "spread_damage", "lc_method", "wl", nullptr};

    if( !PyArg_ParseTupleAndKeywords(
        empty, kwargs, "d|iddi$ppppiiO", (char**)kw,
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
        &wl               // O
    ) )
    {
        Py_DECREF( empty );
        return 0;
    }
    Py_DECREF( empty );

    // Parameters of the SN curve, if defined
    if( wl && wl != Py_None )
    {
        Rainflow::rfc_wl_param_s wl_param = {0};

        if( get_dict_wl( wl, "wl", wl_param, &wl_extended_def ) )
        {
            wl_sx = wl_param.sx;
            wl_nx = wl_param.nx;
            wl_sd = wl_param.sd;
            wl_nd = wl_param.nd;
            wl_k  = wl_param.k;
            wl_k2 = wl_param.k2;
            wl_omission = wl_param.omission;
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

    if( !wl_extended_def )
    {
        if( !rf->wl_init_modified( wl_sx, wl_nx, wl_k, wl_k2 ) )
        {
            PyErr_Format( PyExc_RuntimeError, "Rainflow initialization error (%s)", rfc_err_str( rf->error_get() ) );
            return 0;
        }
    }
    else
    {
        Rainflow::rfc_wl_param_s wl = {0};

        wl.sd = wl_sd;
        wl.nd = wl_nd;
        wl.sx = wl_sx;
        wl.nx = wl_nx;
        wl.k  = wl_k;
        wl.k2 = wl_k2;
        wl.omission = wl_omission;

        if( !rf->wl_init_any( &wl ) )
        {
            PyErr_Format( PyExc_RuntimeError, "Rainflow initialization error (%s)", rfc_err_str( rf->error_get() ) );
            return 0;
        }
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


/**
 * @brief Performs the Rainflow cycle counting algorithm on input data.
 *
 * This function executes the full Rainflow counting process using a configured
 * `Rainflow` instance and a series of numerical input data points. It feeds
 * the data into the Rainflow engine, extracts the raw residuum, and handles
 * a specific post-processing step to remove pending (incomplete) cycles
 * if applicable.
 *
 * @param rf Pointer to a configured `Rainflow` object.
 * @param data Pointer to an array of `npy_double` values (input signal).
 * @param len Number of data points in the input array.
 * @param res_method Residual method to apply when finalizing results.
 * @param residuum_raw Output vector that will contain the raw residuum data.
 *
 * @return `1` on success; `0` if any part of the counting process fails.
 *
 * @throws Raises a Python `RuntimeError` with an explanatory message if feeding,
 *         finalizing, or retrieving results from the `Rainflow` instance fails.
 *
 * @note Post-processing step attempts to detect and remove a specific form
 *       of incomplete cycle pattern at the end of the signal (based on
 *       amplitude symmetry checks).
 */
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
            residuum_raw.erase(residuum_raw.end() - 3, residuum_raw.end() - 1);
        }
    }

    if( !rf->finalize( res_method ) ) goto fail;

    return 1;
fail:
    PyErr_Format( PyExc_RuntimeError, "Error while counting (%s)", rfc_err_str( rf->error_get() ) );
    return 0;
}


/**
 * @brief Prepares and formats the results of the Rainflow counting process as Python objects.
 *
 * This function extracts all computed outputs from a `Rainflow` instance and packages them
 * into a Python dictionary for return to the Python layer. The dictionary includes range pair
 * counts, level crossings, turning points, raw and finalized residue data, rainflow matrix,
 * damage history, and the impaired Wöhler (SN curve) parameters.
 *
 * @param rf Pointer to a `Rainflow` instance containing the counted results.
 * @param data_len Length of the original input data (used for sizing damage history array).
 * @param res_method The method used to handle residuals (used in finalization).
 * @param residuum_raw The raw residue data vector produced during counting.
 * @param ret Output parameter: pointer to a `PyObject*` where the resulting dictionary will be stored.
 *
 * @return `1` on success; `0` if memory allocation or internal retrieval fails.
 *
 * @throws Raises a Python `MemoryError` or `RuntimeError` with descriptive messages
 *         if result preparation fails at any stage.
 *
 * @details The returned Python dictionary includes:
 * - `"damage"`: Total fatigue damage (float)
 * - `"rp"`: Range pairs (NumPy array of shape [class_count, 2])
 * - `"lc"`: Level crossings (NumPy array of shape [class_count, 2])
 * - `"tp"`: Turning points with positions, values, and damage (NumPy array)
 * - `"res_raw"`: Raw residue data (NumPy array)
 * - `"res"`: Finalized residue data (NumPy array)
 * - `"rfm"`: Rainflow matrix (NumPy array of shape [class_count, class_count])
 * - `"dh"`: Damage history per data point (NumPy array of length `data_len`)
 * - `"wl_miner_consistent"`: Impaired Wöhler curve parameters (Python dict)
 */
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
    PyDict_SetItemString( *ret, "wl_miner_consistent", obj );

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


/**
 * @brief Main entry point for performing Rainflow cycle counting from Python.
 *
 * This function receives a Python object containing time series data (typically a NumPy array)
 * and an optional set of keyword arguments that configure the Rainflow counting behavior.
 * It validates and converts the input data, parses all keyword arguments, performs the counting
 * operation, and returns a dictionary of results to Python.
 *
 * @param self Unused. Included for compatibility with the Python C API.
 * @param args Tuple of positional arguments. Expects one argument: the input time series.
 * @param kwargs Dictionary of keyword arguments for configuration (see `parse_rfc_kwargs`).
 *
 * @return A Python dictionary containing Rainflow analysis results (or `nullptr` on failure).
 *
 * @throws Raises Python exceptions (`TypeError`, `ValueError`, `RuntimeError`, etc.) if:
 * - The input data is invalid or cannot be converted to a NumPy array.
 * - The keyword arguments are invalid or inconsistent.
 * - The Rainflow process fails (e.g., memory error or invalid configuration).
 *
 * @details This function delegates to the following components:
 * - `convert_to_numpy_array`: Converts and validates the input array.
 * - `parse_rfc_kwargs`: Parses keyword arguments and initializes the `Rainflow` instance.
 * - `do_rainflow`: Performs the actual counting and residual handling.
 * - `prepare_results`: Formats the output into a Python dictionary.
 *
 * @note Input must be a 1D NumPy-compatible sequence of double precision floats.
 *       Returns `nullptr` and sets a Python exception on failure.
 */
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


/**
 * @brief Computes cumulative fatigue damage based on range pair data and Wöhler (S-N) curve parameters.
 *
 * This function calculates the total fatigue damage using the supplied amplitude values (`Sa`) and
 * corresponding cycle counts (`counts`). Optionally, it uses a dictionary of Wöhler parameters and
 * a selected damage calculation method. The range pairs are sorted by amplitude before computation.
 *
 * @param self Unused. Included for compatibility with the Python C API.
 * @param args Tuple of required arguments:
 *   - `Sa`: A NumPy array of stress amplitudes.
 *   - `counts`: A NumPy array of corresponding cycle counts.
 * @param kwargs Optional keyword arguments:
 *   - `wl`: Dictionary of Wöhler parameters (`sd`, `nd`, `k`, `k2`, `sx`, `nx`, `omission`).
 *   - `method`: Integer code (0–3) for selecting the damage calculation method.
 *
 * @return A Python float representing the computed cumulative damage, or `nullptr` on failure.
 *
 * @throws Raises Python exceptions if:
 * - Input arrays are not convertible to NumPy arrays.
 * - Arrays differ in length or contain invalid (e.g., negative) values.
 * - Wöhler dictionary is malformed or contains non-numeric entries.
 * - An invalid method index is provided.
 *
 * @note Supported damage calculation methods correspond to internal `Rainflow::RFC_RP_DAMAGE_CALC_METHOD_*` enums.
 *       The Wöhler parameters default to `sx=1000`, `nx=1e7`, `k=k2=5` if not provided.
 */
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
        PyErr_SetString(PyExc_ValueError, "`damage_method` must be in range 0 to 3.");
        return nullptr;
    }

    if( Py_wl && Py_wl != Py_None )
    {
        if( !get_dict_wl( Py_wl, "wl", wl ) ) return nullptr;
    }
    else
    {
        wl.k = wl.k2 = 5;
        wl.sx = 1e3;
        wl.nx = 1e7;
        wl.sd = 0.0;
        wl.nd = DBL_MAX;
    }

    if( convert_to_numpy_array(Py_Sa, &arr, &buffer, &size) )
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

    if( convert_to_numpy_array(Py_counts, &arr, &buffer, &size) )
    {
        if( vec_sa.size() != (size_t)size )
        {
            PyErr_SetString( PyExc_ValueError, "`Sa` and `counts` must be of same size." );
        }
        else {
            class_count = (unsigned int)size;
            for( auto i: vec_sorted_indices )
            {
                if( buffer[i] < 0 )
                {
                    PyErr_SetString(PyExc_ValueError, "Negative values in `counts`.");
                    break;
                }
                vec_counts.push_back( (Rainflow::rfc_counts_t)buffer[i] );
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


/**
 * @brief Initializes the `rfcnt` Python module.
 *
 * This is the required entry point for Python’s C extension mechanism. It creates and returns
 * the module object, sets up method definitions, and ensures that the NumPy C API is properly
 * initialized via `import_array()`.
 *
 * @return A pointer to the initialized Python module object, or `nullptr` if initialization fails.
 *
 * @note This function must be named `PyInit_<modulename>` (here, `rfcnt`) to comply with Python 3’s
 *       module initialization requirements.
 *
 * @warning Failing to call `import_array()` will result in crashes or undefined behavior when using
 *          NumPy C functions.
 */
PyMODINIT_FUNC PyInit_rfcnt( void ) {
    PyObject* mod = PyModule_Create( &mymath_module );
    // Initialize numpy
    import_array();
    return mod;
}
