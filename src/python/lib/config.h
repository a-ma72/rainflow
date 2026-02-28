#define ON  1
#define OFF 0

#ifndef RFC_VERSION_MAJOR
  #define RFC_VERSION_MAJOR         "0"
  #define RFC_VERSION_MINOR         "5"
  #define RFC_VERSION_PATCH         "2"
  #define RFC_USE_INTEGRAL_COUNTS    OFF
  #define RFC_USE_HYSTERESIS_FILTER  ON
  #define RFC_MINIMAL                OFF
  #define RFC_TP_SUPPORT             ON
  #define RFC_HCM_SUPPORT            ON
  #define RFC_ASTM_SUPPORT           ON
  #define RFC_USE_DELEGATES          ON
  #define RFC_GLOBAL_EXTREMA         ON
  #define RFC_DAMAGE_FAST            ON
  #define RFC_DH_SUPPORT             ON
  #define RFC_AT_SUPPORT             ON
  #define RFC_AR_SUPPORT             ON
  #define RFC_DEBUG_FLAGS            OFF
  #define RFC_EXPORT_MEX             ON
  #define RFC_EXPORT_PY              ON
  #define RFC_VALUE_TYPE             double
#endif /*RFC_VERSION_MAJOR*/
