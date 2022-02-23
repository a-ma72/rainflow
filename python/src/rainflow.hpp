/*
 *
 *   |                     .-.
 *   |                    /   \
 *   |     .-.===========/     \         .-.
 *   |    /   \         /       \       /   \
 *   |   /     \       /         \     /     \         .-.
 *   +--/-------\-----/-----------\---/-------\-------/---\
 *   | /         \   /             '-'=========\     /     \   /
 *   |/           '-'                           \   /       '-'
 *   |                                           '-'
 *          ____  ___    _____   __________    ____ _       __
 *         / __ \/   |  /  _/ | / / ____/ /   / __ \ |     / /
 *        / /_/ / /| |  / //  |/ / /_  / /   / / / / | /| / /
 *       / _, _/ ___ |_/ // /|  / __/ / /___/ /_/ /| |/ |/ /
 *      /_/ |_/_/  |_/___/_/ |_/_/   /_____/\____/ |__/|__/
 *
 *    Rainflow Counting Algorithm (4-point-method), C99 compliant
 *
 */


/*
    Woehler formula :    (Sa/SD)^-|k| == n/ND
    Basquin formula :    C            == n * Sa^b       (e.g. C = 2e21 if SD=1e3, ND=2e6 and b=5)

    Simplified formula to calculate damage (Miner elementary):
    Fatigue strength:                    SD (=1e3, e.g.)
    Cycles @ SD:                         ND (=2e6, e.g.)
    Woehler slope:                       k  (=5, e.g.)
    Stress amplitude for class i:        Sa_i = ABS(from_i – to_i) * class_width/2
    Cycle counts for class i:            h_i
    Partial damage for class i:          D_i  = h_i/ND * (Sa_i/SD) ^ ABS(k)
                                         D_i  = h_i * Sa ^ b / C
    Damage for entire histogram:         D    = sum( D_i )
*/



#pragma once

#include <vector>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include "rainflow.h"

#pragma pack(push, 1)

/* Suppose correct configuration */
static
char compiler_assert_rfc_config
[
  !RFC_MINIMAL              &&
   RFC_TP_SUPPORT           &&
   RFC_HCM_SUPPORT          &&
   RFC_ASTM_SUPPORT         &&
   RFC_USE_DELEGATES        &&
   RFC_GLOBAL_EXTREMA       &&
   RFC_DAMAGE_FAST          &&
   RFC_DH_SUPPORT           &&
   RFC_AT_SUPPORT           &&
   RFC_AR_SUPPORT
];

#ifndef CALLOC
#define CALLOC calloc
#endif
#ifndef REALLOC
#define REALLOC realloc
#endif
#ifndef FREE
#define FREE free
#endif

namespace RF = RFC_CPP_NAMESPACE;

template< class T = std::vector<RF::rfc_value_tuple_s> > class RainflowT;

template< class T >
static
void * rfc_mem_alloc_default( void *ptr, size_t num, size_t size, int aim )
{
    return RainflowT<T>::mem_alloc( ptr, num, size, (typename RainflowT<T>::rfc_mem_aim_e) aim );
}

#ifndef RFC_MEM_ALLOC
#define RFC_MEM_ALLOC rfc_mem_alloc_default<T>
#endif


/* Support external turning point storage
 *
 * RFC_TP_STORAGE defines the storage container for turning points 
 */
#ifdef RFC_TP_STORAGE

/* C delegates */
extern "C"
{
    static bool  rfc_storage_tp_set           ( RF::rfc_ctx_s* ctx, size_t tp_pos, RF::rfc_value_tuple_s *tp );
    static bool  rfc_storage_tp_get           ( RF::rfc_ctx_s* ctx, size_t tp_pos, RF::rfc_value_tuple_s **tp );
    static bool  rfc_storage_tp_inc_damage    ( RF::rfc_ctx_s *ctx, size_t tp_pos, double damage );
}

typedef RainflowT<RFC_TP_STORAGE> Rainflow;
#else
typedef RainflowT<> Rainflow;

#endif /*RFC_TP_STORAGE*/


template< class T >
class RainflowT
{
public:

    /* Memory allocation aim info */
    enum rfc_mem_aim
    {
        RFC_MEM_AIM_TEMP                        =  RF::RFC_MEM_AIM_TEMP,                        /**< Error on accessing memory for temporary storage */
        RFC_MEM_AIM_RESIDUE                     =  RF::RFC_MEM_AIM_RESIDUE,                     /**< Error on accessing memory for residue */
        RFC_MEM_AIM_MATRIX                      =  RF::RFC_MEM_AIM_MATRIX,                      /**< Error on accessing memory for rf matrix */
        RFC_MEM_AIM_RP                          =  RF::RFC_MEM_AIM_RP,                          /**< Error on accessing memory for range pair counting */
        RFC_MEM_AIM_LC                          =  RF::RFC_MEM_AIM_LC,                          /**< Error on accessing memory for level crossing */
        RFC_MEM_AIM_TP                          =  RF::RFC_MEM_AIM_TP,                          /**< Error on accessing memory for turning points */
        RFC_MEM_AIM_DLUT                        =  RF::RFC_MEM_AIM_DLUT,                        /**< Error on accessing memory for damage look-up table */
        RFC_MEM_AIM_HCM                         =  RF::RFC_MEM_AIM_HCM,                         /**< Error on accessing memory for HCM algorithm */
        RFC_MEM_AIM_DH                          =  RF::RFC_MEM_AIM_DH,                          /**< Error on accessing memory for damage history */
        RFC_MEM_AIM_RFM_ELEMENTS                =  RF::RFC_MEM_AIM_RFM_ELEMENTS,                /**< Error on accessing memory for rf matrix elements */
    };


    /* Flags */
    enum rfc_flags
    {
        RFC_FLAGS_DEFAULT                       = RF::RFC_FLAGS_DEFAULT,                         
        RFC_FLAGS_COUNT_RFM                     = RF::RFC_FLAGS_COUNT_RFM,                      /**< Count into rainflow matrix */
        RFC_FLAGS_COUNT_DAMAGE                  = RF::RFC_FLAGS_COUNT_DAMAGE,                   /**< Count damage */
        RFC_FLAGS_COUNT_DH                      = RF::RFC_FLAGS_COUNT_DH,                       /**< Spread damage */
        RFC_FLAGS_COUNT_RP                      = RF::RFC_FLAGS_COUNT_RP,                       /**< Count into range pair */
        RFC_FLAGS_COUNT_LC_UP                   = RF::RFC_FLAGS_COUNT_LC_UP,                    /**< Count into level crossing (only rising slopes) */
        RFC_FLAGS_COUNT_LC_DN                   = RF::RFC_FLAGS_COUNT_LC_DN,                    /**< Count into level crossing (only falling slopes) */
        RFC_FLAGS_COUNT_LC                      = RF::RFC_FLAGS_COUNT_LC,                       /**< Count into level crossing (all slopes) */
        RFC_FLAGS_COUNT_MK                      = RF::RFC_FLAGS_COUNT_MK,                       /**< Live damage counter (Miner consequent) */
        RFC_FLAGS_ENFORCE_MARGIN                = RF::RFC_FLAGS_ENFORCE_MARGIN,                 /**< Enforce first and last data point are turning points */
        RFC_FLAGS_COUNT_ALL                     = RF::RFC_FLAGS_COUNT_ALL,                      /**< Count all */
        RFC_FLAGS_TPPRUNE_PRESERVE_POS          = RF::RFC_FLAGS_TPPRUNE_PRESERVE_POS,           /**< Preserve stream position information on pruning */
        RFC_FLAGS_TPPRUNE_PRESERVE_RES          = RF::RFC_FLAGS_TPPRUNE_PRESERVE_RES,           /**< Preserve turning points that exist in resiude on pruning */
        RFC_FLAGS_TPAUTOPRUNE                   = RF::RFC_FLAGS_TPAUTOPRUNE,                    /**< Automatic prune on tp */
        RFC_FLAGS_AUTORESIZE                    = RF::RFC_FLAGS_AUTORESIZE,                     /**< Automatically resize buffers for rp, lc, and rfm */
    };


    enum rfc_debug_flags
    {
        RFC_FLAGS_LOG_CLOSED_CYCLES             = RF::RFC_FLAGS_LOG_CLOSED_CYCLES,              /**< Log closed cycles */
    };



    /* See RFC_damage_from_rp() */
    enum rfc_rp_damage_method
    {
        RFC_RP_DAMAGE_CALC_METHOD_DEFAULT       = RF::RFC_RP_DAMAGE_CALC_METHOD_DEFAULT,        /**< Use Woehler parameters from rfc_ctx */
        RFC_RP_DAMAGE_CALC_METHOD_ELEMENTAR     = RF::RFC_RP_DAMAGE_CALC_METHOD_ELEMENTAR,      /**< Use Woehler parameters from rfc_ctx, but as Miner elementar */
        RFC_RP_DAMAGE_CALC_METHOD_MODIFIED      = RF::RFC_RP_DAMAGE_CALC_METHOD_MODIFIED,       /**< Use Woehler parameters from rfc_ctx, but as Miner modified */
        RFC_RP_DAMAGE_CALC_METHOD_CONSEQUENT    = RF::RFC_RP_DAMAGE_CALC_METHOD_CONSEQUENT,     /**< Use Woehler parameters from rfc_ctx, but as Miner consequent */
    };


    enum rfc_state
    {
        RFC_STATE_INIT0                         = RF::RFC_STATE_INIT0,                          /**< Initialized with zeros */
        RFC_STATE_INIT                          = RF::RFC_STATE_INIT,                           /**< Initialized, memory allocated */
        RFC_STATE_BUSY                          = RF::RFC_STATE_BUSY,                           /**< In counting state */
        RFC_STATE_BUSY_INTERIM                  = RF::RFC_STATE_BUSY_INTERIM,                   /**< In counting state, having still one interim turning point (not included) */
        RFC_STATE_FINALIZE                      = RF::RFC_STATE_FINALIZE,                       /**< Finalizing */
        RFC_STATE_FINISHED                      = RF::RFC_STATE_FINISHED,                       /**< Counting finished, memory still allocated */
        RFC_STATE_ERROR                         = RF::RFC_STATE_ERROR,                          /**< An error occurred */
    };


    enum rfc_error
    {
        RFC_ERROR_UNEXP                         = RF::RFC_ERROR_UNEXP,                          /**< Unexpected error */
        RFC_ERROR_NOERROR                       = RF::RFC_ERROR_NOERROR,                        /**< No error */
        RFC_ERROR_INVARG                        = RF::RFC_ERROR_INVARG,                         /**< Invalid arguments passed */
        RFC_ERROR_UNSUPPORTED                   = RF::RFC_ERROR_UNSUPPORTED,                    /**< Unsupported feature */
        RFC_ERROR_MEMORY                        = RF::RFC_ERROR_MEMORY,                         /**< Error on memory allocation */
        RFC_ERROR_TP                            = RF::RFC_ERROR_TP,                             /**< Error while processing turning points */
        RFC_ERROR_AT                            = RF::RFC_ERROR_AT,                             /**< Error while amplitude transformation */
        RFC_ERROR_DH_BAD_STREAM                 = RF::RFC_ERROR_DH_BAD_STREAM,                  /**< Input stream must be unique */
        RFC_ERROR_DH                            = RF::RFC_ERROR_DH,                             /**< Error while damage history calculation/access */
        RFC_ERROR_LUT                           = RF::RFC_ERROR_LUT,                            /**< Error while accessing look up tables */
        RFC_ERROR_DATA_OUT_OF_RANGE             = RF::RFC_ERROR_DATA_OUT_OF_RANGE,              /**< Input data leaves classrange */
        RFC_ERROR_DATA_INCONSISTENT             = RF::RFC_ERROR_DATA_INCONSISTENT,              /**< Processed data is inconsistent (internal error) */
    };


    enum rfc_counting_method
    {
        RFC_COUNTING_METHOD_DELEGATED           = RF::RFC_COUNTING_METHOD_DELEGATED,            /**< Method must be implemented via delegator, see member cycle_find_fcn */
        RFC_COUNTING_METHOD_NONE                = RF::RFC_COUNTING_METHOD_NONE,                 /**< No counting */
        RFC_COUNTING_METHOD_4PTM                = RF::RFC_COUNTING_METHOD_4PTM,                 /**< 4 point algorithm (default) */
        RFC_COUNTING_METHOD_HCM                 = RF::RFC_COUNTING_METHOD_HCM,                  /**< 3 point algorithm, Clormann/Seeger (HCM) method */
        RFC_COUNTING_METHOD_ASTM                = RF::RFC_COUNTING_METHOD_ASTM,                 /**< 3 point algorithm, ASTM Standard E 1049 */
        RFC_COUNTING_METHOD_COUNT               = RF::RFC_COUNTING_METHOD_COUNT,                /**< Number of options */
    };


    enum rfc_res_method
    {
        /* Don't change order! */
        RFC_RES_NONE                            = RF::RFC_RES_NONE,                             /**< No residual method */
        RFC_RES_IGNORE                          = RF::RFC_RES_IGNORE,                           /**< Ignore residue (same as RFC_RES_NONE) */
        RFC_RES_NO_FINALIZE                     = RF::RFC_RES_NO_FINALIZE,                      /**< Don't finalize data stream */
        RFC_RES_DISCARD                         = RF::RFC_RES_DISCARD,                          /**< Discard residue (empty residue) */
        RFC_RES_HALFCYCLES                      = RF::RFC_RES_HALFCYCLES,                       /**< ASTM */
        RFC_RES_FULLCYCLES                      = RF::RFC_RES_FULLCYCLES,                       /**< Count half cycles as full cycles */
        RFC_RES_CLORMANN_SEEGER                 = RF::RFC_RES_CLORMANN_SEEGER,                  /**< Clormann/Seeger method */
        RFC_RES_REPEATED                        = RF::RFC_RES_REPEATED,                         /**< Repeat residue and count closed cycles */
        RFC_RES_RP_DIN45667                     = RF::RFC_RES_RP_DIN45667,                      /**< Count residue according to range pair in DIN-45667 */
        RFC_RES_COUNT                           = RF::RFC_RES_COUNT,                            /**< Number of options */
    };


    enum rfc_sd_method
    {
        RFC_SD_NONE                             = RF::RFC_SD_NONE,                              /**< No spread damage calculation */
        RFC_SD_HALF_23                          = RF::RFC_SD_HALF_23,                           /**< Equally split damage between P2 and P3 */
        RFC_SD_RAMP_AMPLITUDE_23                = RF::RFC_SD_RAMP_AMPLITUDE_23,                 /**< Spread damage according to amplitude over points between P2 and P3 */
        RFC_SD_RAMP_DAMAGE_23                   = RF::RFC_SD_RAMP_DAMAGE_23,                    /**< Spread damage linearly over points between P2 and P3 */
        RFC_SD_RAMP_AMPLITUDE_24                = RF::RFC_SD_RAMP_AMPLITUDE_24,                 /**< Spread damage according to amplitude over points between P2 and P4 */  
        RFC_SD_RAMP_DAMAGE_24                   = RF::RFC_SD_RAMP_DAMAGE_24,                    /**< Spread damage linearly over points between P2 and P4 */
        RFC_SD_FULL_P2                          = RF::RFC_SD_FULL_P2,                           /**< Assign damage to P2 */
        RFC_SD_FULL_P3                          = RF::RFC_SD_FULL_P3,                           /**< Assign damage to P3 */
        RFC_SD_TRANSIENT_23                     = RF::RFC_SD_TRANSIENT_23,                      /**< Spread damage transient according to amplitude over points between P2 and P3 */
        RFC_SD_TRANSIENT_23c                    = RF::RFC_SD_TRANSIENT_23c,                     /**< Spread damage transient according to amplitude over points between P2 and P4 only until cycle is closed */
        RFC_SD_COUNT                            = RF::RFC_SD_COUNT,                             /**< Number of options */
    };


    enum rfc_lc_count_method
    {
        RFC_LC_COUNT_METHOD_SLOPES_UP           = RF::RFC_LC_COUNT_METHOD_SLOPES_UP,            /**< Count on rising slopes only (default) */
        RFC_LC_COUNT_METHOD_SLOPES_DOWN         = RF::RFC_LC_COUNT_METHOD_SLOPES_DOWN,          /**< Count on falling slopes only */
        RFC_LC_COUNT_METHOD_SLOPES_ALL          = RF::RFC_LC_COUNT_METHOD_SLOPES_ALL,           /**< Count on rising AND falling slopes */
    };


    /* Typedefs */
    typedef                 RF::rfc_value_t         rfc_value_t;                                /** Input data value type */
    typedef                 RF::rfc_counts_t        rfc_counts_t;                               /** Type of counting values */
    typedef                 RF::rfc_value_tuple_s   rfc_value_tuple_s;                          /** Tuple of value and index position */
    typedef                 RF::rfc_ctx_s           rfc_ctx_s;                                  /** Forward declaration (rainflow context) */
    typedef                 RF::rfc_class_param     rfc_class_param_s;                          /** Class parameters (width, offset, count) */
    typedef                 RF::rfc_wl_param        rfc_wl_param_s;                             /** Woehler curve parameters (sd, nd, k, k2, omission) */
    typedef                 RF::rfc_rfm_item        rfc_rfm_item_s;                             /** Rainflow matrix element */
    typedef     enum        rfc_mem_aim             rfc_mem_aim_e;                              /** Memory accessing mode */
    typedef     enum        rfc_flags               rfc_flags_e;                                /** Flags, see RFC_FLAGS... */
    typedef     enum        rfc_state               rfc_state_e;                                /** Counting state, see RFC_STATE... */
    typedef     enum        rfc_error               rfc_error_e;                                /** Recent error, see RFC_ERROR... */
    typedef     enum        rfc_res_method          rfc_res_method_e;                           /** Method when count residue into matrix, see RFC_RES... */
    typedef     enum        rfc_counting_method     rfc_counting_method_e;                      /** Counting method, see RFC_COUNTING... */
    typedef     enum        rfc_rp_damage_method    rfc_rp_damage_method_e;                     /** Method when calculating damage from range pair counting, see RFC_RP_DAMAGE_CALC_METHOD... */
    typedef     enum        rfc_sd_method           rfc_sd_method_e;                            /** Spread damage method, see RFC_SD... */
    typedef     enum        rfc_lc_count_method     rfc_lc_count_method_e;                      /** Controls which slopes to take into account, when doing the level crossing counting */

    typedef     std::vector<double>                 rfc_double_v;                               /** Vector of double */
    typedef     std::vector<rfc_value_t>            rfc_value_v;                                /** Vector of values */
    typedef     std::vector<rfc_counts_t>           rfc_counts_v;                               /** Vector of counts */
    typedef     std::vector<rfc_rfm_item_s>         rfc_rfm_item_v;                             /** Vector of rainflow matrix items */

    typedef     T                                   rfc_tp_storage;                             /** Rainflow turning points storage */


    /* Memory allocation functions typedef */
    typedef     void *   ( *rfc_mem_alloc_fcn_t )   ( void *, size_t num, size_t size, rfc_mem_aim_e aim );     /** Memory allocation functor */

    /* Core function wrapper */
    bool            init                    ( unsigned class_count, rfc_value_t class_width, rfc_value_t class_offset, 
                                              rfc_value_t hysteresis, rfc_flags_e flags = RFC_FLAGS_DEFAULT );
    rfc_state_e     state_get               () const { return (rfc_state_e)RF::RFC_state_get( &ctx_get() ); }
    rfc_error_e     error_get               () const { return (rfc_error_e)RF::RFC_error_get( &ctx_get() ); }
    bool            wl_init_elementary      ( double sx, double nx, double k );
    bool            wl_init_original        ( double sd, double nd, double k );
    bool            wl_init_modified        ( double sx, double nx, double k, double k2 );
    bool            wl_init_any             ( const rfc_wl_param_s* );
    bool            clear_counts            ();
    bool            deinit                  ();
    bool            feed                    ( const rfc_value_t* data, size_t count );
    bool            cycle_process_counts    ( rfc_value_t from_val, rfc_value_t to_val, rfc_flags_e flags );
    bool            feed_scaled             ( const rfc_value_t* data, size_t count, double factor );
    bool            feed_tuple              ( rfc_value_tuple_s *data, size_t count );
    bool            finalize                ( rfc_res_method_e residual_method = RFC_RES_IGNORE );
    /* Functions on rainflow matrix */           
    bool            rfm_make_symmetric      ();
    bool            rfm_non_zeros           ( unsigned *count ) const;
    bool            rfm_get                 ( rfc_rfm_item_s **buffer, unsigned *count ) const;
    bool            rfm_set                 ( const rfc_rfm_item_s *buffer, unsigned count, bool add_only );
    bool            rfm_peek                ( rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t *count ) const;
    bool            rfm_poke                ( rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t count, bool add_only );
    bool            rfm_sum                 ( unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, rfc_counts_t *count ) const;
    bool            rfm_damage              ( unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, double *damage ) const;
    bool            rfm_check               () const;
    bool            rfm_refeed              ( rfc_value_t new_hysteresis, const rfc_class_param_s *new_class_param );
/* Functions on histograms */           
    bool            lc_get                  ( rfc_counts_t *lc, rfc_value_t *level ) const;
    bool            lc_from_rfm             ( rfc_counts_t *lc, rfc_value_t *level, const rfc_counts_t *rfm, rfc_flags_e flags ) const;
    bool            lc_from_residue         ( rfc_counts_t *lc, rfc_value_t *level, rfc_flags_e flags ) const;
    bool            rp_get                  ( rfc_counts_t *rp, rfc_value_t *Sa ) const;
    bool            rp_from_rfm             ( rfc_counts_t *rp, rfc_value_t *Sa, const rfc_counts_t *rfm ) const;
    bool            damage                  ( rfc_value_t *damage = NULL, rfc_value_t *damage_residue = NULL ) const;
    bool            damage_from_rp          ( const rfc_counts_t *counts, const rfc_value_t *Sa, double *damage, rfc_rp_damage_method_e rp_calc_type ) const;
    bool            damage_from_rfm         ( const rfc_counts_t *rfm, double *damage ) const;
    /* Woehler curve */
    bool            wl_calc_sx              ( double s0, double n0, double k, double *sx, double nx, double  k2, double  sd, double nd ) const;
    bool            wl_calc_sd              ( double s0, double n0, double k, double  sx, double nx, double  k2, double *sd, double nd ) const;
    bool            wl_calc_k2              ( double s0, double n0, double k, double  sx, double nx, double *k2, double  sd, double nd ) const;
    bool            wl_calc_sa              ( double s0, double n0, double k, double  n,  double *sa ) const;
    bool            wl_calc_n               ( double s0, double n0, double k, double  sa, double *n ) const;
    /* Turning points */
    bool            tp_init_autoprune       ( bool autoprune, size_t size, size_t threshold );
    bool            tp_prune                ( size_t count, rfc_flags_e flags );
    bool            tp_refeed               ( rfc_value_t new_hysteresis, const rfc_class_param_s *new_class_param );
    bool            tp_clear                ();
    /* Residuum */
    bool            res_get                 ( const rfc_value_tuple_s **residue, unsigned *count ) const;
    /* Damage history */
    bool            dh_init                 ( rfc_sd_method_e method, double *dh, size_t dh_cap, bool is_static );
    bool            dh_get                  ( const double **dh, size_t *count ) const;
    /* Amplitude transformation*/
    bool            at_init                 ( const double *Sa, const double *Sm, unsigned count, 
                                              double M, double Sm_rig, double R_rig, bool R_pinned, bool symmetric );
    bool            at_transform            ( double Sa, double Sm, double *Sa_transformed ) const;
    /* Flags */
    bool            flags_set               ( int flags, bool debugging = false, bool overwrite = true );
    bool            flags_unset             ( int flags, bool debugging = false );
    bool            flags_get               ( int *flags, bool debugging = false ) const;
    inline
    rfc_counts_t    full_inc                () const { return m_ctx.full_inc; }
    inline
    rfc_counts_t    half_inc                () const { return m_ctx.half_inc; }
    bool            cls_number              ( rfc_value_t value, unsigned *class_number ) const;
    bool            cls_upper               ( unsigned class_number, rfc_value_t *class_upper ) const;
    bool            cls_mean                ( unsigned class_number, rfc_value_t *class_mean ) const;
    bool            class_count             ( unsigned *class_count ) const;
    bool            class_offset            ( rfc_value_t *class_offset ) const;
    bool            class_width             ( rfc_value_t *class_width ) const;
    bool            hysteresis              ( rfc_value_t *hysteresis ) const;

    /* more C++ specific extensions */
    bool            feed                    ( const std::vector<rfc_value_t> data );
    bool            feed_scaled             ( const std::vector<rfc_value_t> data, double factor );
    bool            rfm_get                 ( rfc_rfm_item_v &buffer ) const;
    bool            rfm_set                 ( const rfc_rfm_item_v &buffer, bool add_only );
    bool            lc_get                  ( rfc_counts_v &lc, rfc_value_v &level ) const;
    bool            lc_from_rfm             ( rfc_counts_v &lc, rfc_value_v &level, const rfc_counts_t *rfm, rfc_flags_e flags ) const;
    bool            lc_from_residue         ( rfc_counts_v &lc, rfc_value_v &level, rfc_flags_e flags ) const;
    bool            rp_get                  ( rfc_counts_v &rp, rfc_value_v &Sa ) const;
    bool            rp_from_rfm             ( rfc_counts_v &rp, rfc_value_v &Sa, const rfc_counts_t *rfm ) const;
    bool            damage_from_rp          ( const rfc_counts_v &counts, const rfc_value_v &Sa, double &damage, rfc_rp_damage_method_e rp_calc_type ) const;
    bool            at_init                 ( const rfc_double_v &Sa, const rfc_double_v &Sm, 
                                              double M, double Sm_rig, double R_rig, bool R_pinned, bool symmetric );
    bool            at_init                 ( double M, double Sm_rig, double R_rig, bool R_pinned );
    bool            at_transform            ( double Sa, double Sm, double &Sa_transformed ) const;
    bool            wl_param_get            ( rfc_wl_param_s &wl_param ) const;

    /* TP storage access */
    inline const
    rfc_tp_storage& tp_storage              () const { return m_tp; }
    inline
    rfc_tp_storage& tp_storage              ()       { return m_tp; }

    /* matrix access */
    inline const
    rfc_counts_t*   rfm_storage             () const { return m_ctx.rfm; }
    inline
    rfc_counts_t*   rfm_storage             ()       { return m_ctx.rfm; }

    /* Delegates */
    bool            tp_set                  ( size_t tp_pos, rfc_value_tuple_s *tp );
    bool            tp_get                  ( size_t tp_pos, rfc_value_tuple_s **tp );
    bool            tp_inc_damage           ( size_t tp_pos, double damage );

    
    /* dtor */     ~RainflowT               () { deinit(); }
    /* ctor */      RainflowT               ()                                                   // Std ctor
    { 
        rfc_ctx_s nil = { sizeof( rfc_ctx_s ) };

        m_ctx = nil;
        m_ctx.mem_alloc = RFC_MEM_ALLOC;  // wrapper calls class method mem_alloc per default
    } 
    /* ctor */      RainflowT               ( rfc_ctx_s&& other ) { ctx_assign( other ); }   // Move ctor
    RainflowT&      operator=               ( rfc_ctx_s&& other ) { ctx_assign( other ); }   // Move assignment

    /* ctx access */
    const
    rfc_ctx_s&      ctx_get                 () const { return m_ctx; }
    rfc_ctx_s&      ctx_get                 ()       { return m_ctx; }
    void            ctx_assign              ( rfc_ctx_s& ctx )
    { 
        if( ctx.internal.obj != this ) 
        { 
            rfc_ctx_s nil = { sizeof( RF::rfc_ctx_s ) };

            (void)deinit(); 
            m_ctx = ctx;
            m_ctx.internal.obj = this;  // Take ownership and custody
            ctx = nil; 
        } 
    }
    void            ctx_assign              ( rfc_ctx_s&& ctx )  // Move assignment
    { 
        if( ctx.internal.obj != this ) 
        { 
            (void)deinit(); 
            m_ctx = ctx;
            m_ctx.internal.obj = this;  // Take ownership and custody
        } 
    }

    /* Memory allocator */
    static
    void*           mem_alloc               ( void *ptr, size_t num, size_t size, rfc_mem_aim_e aim );

private:
    void            ctx_assign              ( const rfc_ctx_s& );   // Inhibit assign on const ctx
                    RainflowT               ( const rfc_ctx_s& );   // Inhibit copy ctor on const ctx
                    RainflowT               ( const RainflowT& );       // Inhibit copy ctor on (non-)const RainflowT
    RainflowT&      operator=               ( const rfc_ctx_s& );       // Inhibit copy assignment on const ctx
    RainflowT&      operator=               ( const RainflowT& );       // Inhibit copy assignment on (non-)const RainflowT

protected:
    rfc_ctx_s       m_ctx;
    rfc_tp_storage  m_tp;
};




template< class T >
bool RainflowT<T>::init( unsigned class_count, rfc_value_t class_width, rfc_value_t class_offset, 
                         rfc_value_t hysteresis, rfc_flags_e flags )
{
    bool ok;

    ok = RF::RFC_init( &m_ctx, class_count, class_width, class_offset, hysteresis, (RF::rfc_flags_e)flags );

    if( ok )
    {
        m_ctx.internal.obj          = this;
#ifdef RFC_TP_STORAGE
        m_ctx.tp_set_fcn            = rfc_storage_tp_set;
        m_ctx.tp_get_fcn            = rfc_storage_tp_get;
        m_ctx.tp_inc_damage_fcn     = rfc_storage_tp_inc_damage;
#endif /*RFC_TP_STORAGE*/
    }

    return ok;
}


template< class T >
bool RainflowT<T>::wl_init_elementary( double sx, double nx, double k )
{
    return RF::RFC_wl_init_elementary( &m_ctx, sx, nx, k );
}


template< class T >
bool RainflowT<T>::wl_init_original( double sd, double nd, double k )
{
    return RF::RFC_wl_init_original( &m_ctx, sd, nd, k );
}


template< class T >
bool RainflowT<T>::wl_init_modified( double sx, double nx, double k, double k2 )
{
    return RF::RFC_wl_init_modified( &m_ctx, sx, nx, k, k2 );
}


template< class T >
bool RainflowT<T>::wl_init_any( const rfc_wl_param_s* wl_param )
{
    return RF::RFC_wl_init_any( &m_ctx, (const RF::rfc_wl_param_s*) wl_param );
}


template< class T >
bool RainflowT<T>::clear_counts()
{
    return RF::RFC_clear_counts( &m_ctx );
}


template< class T >
bool RainflowT<T>::deinit()
{
    if( !m_ctx.internal.obj )
    {
        return true;
    }

    if( m_ctx.internal.obj != this )
    {
        // Don't have ownership
        return false;
    }
    else
    {
        return RF::RFC_deinit( &m_ctx );
    }
}


template< class T >
bool RainflowT<T>::feed( const rfc_value_t* data, size_t count )
{
    return RF::RFC_feed( &m_ctx, (const RF::rfc_value_t*)data, count );
}


template< class T >
bool RainflowT<T>::cycle_process_counts( rfc_value_t from_val, rfc_value_t to_val, rfc_flags_e flags )
{
    return RF::RFC_cycle_process_counts( &m_ctx, (RF::rfc_value_t)from_val, (RF::rfc_value_t)to_val, (RF::rfc_flags_e)flags );
}


template< class T >
bool RainflowT<T>::feed_scaled( const rfc_value_t* data, size_t count, double factor )
{
    return RF::RFC_feed_scaled( &m_ctx, (const RF::rfc_value_t*)data, count, factor );
}


template< class T >
bool RainflowT<T>::feed_tuple( rfc_value_tuple_s *data, size_t count )
{
    return RF::RFC_feed_tuple( &m_ctx, (RF::rfc_value_tuple_s *)data, count );
}


template< class T >
bool RainflowT<T>::finalize( rfc_res_method_e residual_method )
{
    return RF::RFC_finalize( &m_ctx, (RF::rfc_res_method_e)residual_method );
}


template< class T >
bool RainflowT<T>::rfm_make_symmetric()
{
    return RF::RFC_rfm_make_symmetric( &m_ctx );
}


template< class T >
bool RainflowT<T>::rfm_non_zeros( unsigned *count ) const
{
    return RF::RFC_rfm_non_zeros( &m_ctx, count );
}


template< class T >
bool RainflowT<T>::rfm_get( rfc_rfm_item_s **buffer, unsigned *count ) const
{
    return RF::RFC_rfm_get( &m_ctx, (RF::rfc_rfm_item_s **)buffer, count );
}


template< class T >
bool RainflowT<T>::rfm_set( const rfc_rfm_item_s *buffer, unsigned count, bool add_only )
{
    return RF::RFC_rfm_set( &m_ctx, (const RF::rfc_rfm_item_s *)buffer, count, add_only );
}


template< class T >
bool RainflowT<T>::rfm_peek( rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t *counts ) const
{
    return RF::RFC_rfm_peek( &m_ctx, (RF::rfc_value_t)from_val, (RF::rfc_value_t)to_val, (RF::rfc_counts_t *)counts );
}


template< class T >
bool RainflowT<T>::rfm_poke( rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t counts, bool add_only )
{
    return RF::RFC_rfm_poke( &m_ctx, (RF::rfc_value_t)from_val, (RF::rfc_value_t)to_val, (RF::rfc_counts_t)counts, add_only );
}


template< class T >
bool RainflowT<T>::rfm_sum( unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, rfc_counts_t *count ) const
{
    return RF::RFC_rfm_sum( &m_ctx, from_first, from_last, to_first, to_last, (RF::rfc_counts_t *)count );
}


template< class T >
bool RainflowT<T>::rfm_damage( unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, double *damage ) const
{
    return RF::RFC_rfm_damage( &m_ctx, from_first, from_last, to_first, to_last, damage );
}


template< class T >
bool RainflowT<T>::rfm_refeed( rfc_value_t new_hysteresis, const rfc_class_param_s *new_class_param )
{
    return RF::RFC_rfm_refeed( &m_ctx, new_hysteresis, new_class_param );
}


template< class T >
bool RainflowT<T>::rfm_check() const
{
    return RF::RFC_rfm_check( &m_ctx );
}


template< class T >
bool RainflowT<T>::lc_get( rfc_counts_t *lc, rfc_value_t *level ) const
{
    return RF::RFC_lc_get( &m_ctx, (RF::rfc_counts_t *)lc, (RF::rfc_value_t *)level );
}


template< class T >
bool RainflowT<T>::lc_from_rfm( rfc_counts_t *lc, rfc_value_t *level, const rfc_counts_t *rfm, rfc_flags_e flags ) const
{
    return RF::RFC_lc_from_rfm( &m_ctx, (RF::rfc_counts_t *)lc, (RF::rfc_value_t *)level, (const RF::rfc_counts_t *)rfm, (RF::rfc_flags_e)flags );
}


template< class T >
bool RainflowT<T>::lc_from_residue( rfc_counts_t *lc, rfc_value_t *level, rfc_flags_e flags ) const
{
    return RF::RFC_lc_from_residue( &m_ctx, (RF::rfc_counts_t *)lc, (RF::rfc_value_t *)level, (RF::rfc_flags_e) flags );
}


template< class T >
bool RainflowT<T>::rp_get( rfc_counts_t *rp, rfc_value_t *Sa ) const
{
    return RF::RFC_rp_get( &m_ctx, (RF::rfc_counts_t *)rp, (RF::rfc_value_t *)Sa );
}


template< class T >
bool RainflowT<T>::rp_from_rfm( rfc_counts_t *rp, rfc_value_t *Sa, const rfc_counts_t *rfm ) const
{
    return RF::RFC_rp_from_rfm( &m_ctx, (RF::rfc_counts_t *)rp, (RF::rfc_value_t *)Sa, (const RF::rfc_counts_t *)rfm );
}

template< class T >
bool RainflowT<T>::damage( rfc_value_t *damage, rfc_value_t *damage_residue ) const
{
    return RF::RFC_damage( &m_ctx, damage, damage_residue );
}


template< class T >
bool RainflowT<T>::damage_from_rp( const rfc_counts_t *counts, const rfc_value_t *Sa, double *damage, rfc_rp_damage_method_e rp_calc_type ) const
{
    return RF::RFC_damage_from_rp( &m_ctx, (const RF::rfc_counts_t *)counts, (const RF::rfc_value_t *)Sa, damage, (RF::rfc_rp_damage_method_e)rp_calc_type );
}


template< class T >
bool RainflowT<T>::damage_from_rfm( const rfc_counts_t *rfm, double *damage ) const
{
    return RF::RFC_damage_from_rfm( &m_ctx, (const RF::rfc_counts_t *)rfm, damage );
}


template< class T >
bool RainflowT<T>::wl_calc_sx( double s0, double n0, double k, double *sx, double nx, double  k2, double  sd, double nd ) const
{
    return RF::RFC_wl_calc_sx( &m_ctx, s0, n0, k, sx, nx, k2, sd, nd );
}


template< class T >
bool RainflowT<T>::wl_calc_sd( double s0, double n0, double k, double  sx, double nx, double  k2, double *sd, double nd ) const
{
    return RF::RFC_wl_calc_sd( &m_ctx, s0, n0, k, sx, nx, k2, sd, nd );
}


template< class T >
bool RainflowT<T>::wl_calc_k2( double s0, double n0, double k, double  sx, double nx, double *k2, double  sd, double nd ) const
{
    return RF::RFC_wl_calc_k2( &m_ctx, s0, n0, k, sx, nx, k2, sd, nd );
}


template< class T >
bool RainflowT<T>::wl_calc_sa( double s0, double n0, double k, double  n,  double *sa ) const
{
    return RF::RFC_wl_calc_sa( &m_ctx, s0, n0, k, n,  sa );
}


template< class T >
bool RainflowT<T>::wl_calc_n( double s0, double n0, double k, double  sa, double *n ) const
{
    return RF::RFC_wl_calc_n( &m_ctx, s0, n0, k, sa, n );
}


template< class T >
bool RainflowT<T>::tp_init_autoprune( bool autoprune, size_t size, size_t threshold )
{
    return RF::RFC_tp_init_autoprune( &m_ctx, autoprune, size, threshold );
}


template< class T >
bool RainflowT<T>::tp_prune( size_t count, rfc_flags_e flags )
{
    bool ok;
    
    ok = RF::RFC_tp_prune( &m_ctx, count, (RF::rfc_flags_e) flags );

    m_tp.resize( m_ctx.tp_cnt );

    return ok;
}


template< class T >
bool RainflowT<T>::tp_refeed( rfc_value_t new_hysteresis, const rfc_class_param_s *new_class_param )
{
    bool ok;

    ok = RF::RFC_tp_refeed( &m_ctx, new_hysteresis, new_class_param );

    m_tp.resize( m_ctx.tp_cnt );

    return ok;
}


template< class T >
bool RainflowT<T>::tp_clear()
{
    return RF::RFC_tp_clear( &m_ctx );
}


template< class T >
bool RainflowT<T>::res_get( const rfc_value_tuple_s **residue, unsigned *count ) const
{
    return RF::RFC_res_get( &m_ctx, residue, count );
}


template< class T >
bool RainflowT<T>::dh_init( rfc_sd_method_e method, double *dh, size_t dh_cap, bool is_static )
{
    return RF::RFC_dh_init( &m_ctx, (RF::rfc_sd_method_e)method, dh, dh_cap, is_static );
}


template< class T >
bool RainflowT<T>::dh_get( const double **dh, size_t *count ) const
{
    return RF::RFC_dh_get( &m_ctx, dh, count );
}


template< class T >
bool RainflowT<T>::at_init( const double *Sa, const double *Sm, unsigned count,
                                  double M, double Sm_rig, double R_rig, bool R_pinned, bool symmetric )
{
    return RF::RFC_at_init( &m_ctx, Sa, Sm, count, M, Sm_rig, R_rig, R_pinned, symmetric );
}


template< class T >
bool RainflowT<T>::at_transform( double Sa, double Sm, double *Sa_transformed ) const
{
    return RF::RFC_at_transform( &m_ctx, Sa, Sm, Sa_transformed );
}


template< class T >
bool RainflowT<T>::flags_set( int flags, bool debugging, bool overwrite )
{
    return RF::RFC_flags_set( &m_ctx, flags, debugging, overwrite );
}


template< class T >
bool RainflowT<T>::flags_unset( int flags, bool debugging )
{
    return RF::RFC_flags_unset( &m_ctx, flags, debugging );
}


template< class T >
bool RainflowT<T>::flags_get( int *flags, bool debugging ) const
{
    return RF::RFC_flags_get( &m_ctx, flags, debugging );
}


template< class T >
bool RainflowT<T>::cls_number( rfc_value_t value, unsigned *class_number ) const
{
    return RF::RFC_class_number( &m_ctx, value, class_number );
}


template< class T >
bool RainflowT<T>::cls_upper( unsigned class_number, rfc_value_t *class_upper ) const
{
    return RF::RFC_class_upper( &m_ctx, class_number, class_upper );
}


template< class T >
bool RainflowT<T>::cls_mean( unsigned class_number, rfc_value_t *class_mean ) const
{
    return RF::RFC_class_mean( &m_ctx, class_number, class_mean );
}


template< class T >
bool RainflowT<T>::class_count( unsigned *class_count ) const
{
    return RF::RFC_class_count( &m_ctx, class_count );
}


template< class T >
bool RainflowT<T>::class_offset( rfc_value_t *class_offset ) const
{
    return RF::RFC_class_offset( &m_ctx, class_offset );
}


template< class T >
bool RainflowT<T>::class_width( rfc_value_t *class_width ) const
{
    return RF::RFC_class_width( &m_ctx, class_width );
}


template< class T >
bool RainflowT<T>::hysteresis( rfc_value_t *hysteresis ) const
{
    return RF::RFC_hysteresis( &m_ctx, hysteresis );
}


/* CPP specific extensions */
template< class T >
bool RainflowT<T>::feed( const std::vector<rfc_value_t> data )
{
    return feed( &data[0], data.size() );
}


template< class T >
bool RainflowT<T>::feed_scaled( const std::vector<rfc_value_t> data, double factor )
{
    return feed_scaled( &data[0], data.size(), factor );
}


template< class T >
bool RainflowT<T>::rfm_get( rfc_rfm_item_v &buffer ) const
{
    rfc_rfm_item_s *buffer_ = NULL;
    unsigned        count   = 0;
    bool            ok;

    if( rfm_get( &buffer_, &count ) )
    {
        buffer = rfc_rfm_item_v( buffer_, buffer_ + count );
        ok     = true;
    }
    else
    { 
        ok = false;
    }

    (void)mem_alloc( buffer_, 0, 0, RFC_MEM_AIM_RFM_ELEMENTS );

    return ok;
}


template< class T >
bool RainflowT<T>::rfm_set( const rfc_rfm_item_v &buffer, bool add_only )
{
    return rfm_set( &buffer[0], (unsigned)buffer.size(), add_only );
}


template< class T >
bool RainflowT<T>::lc_get( rfc_counts_v &lc, rfc_value_v &level ) const
{
    lc.resize( m_ctx.class_count );
    level.resize( m_ctx.class_count );

    return lc_get( &lc[0], &level[0] );
}


template< class T >
bool RainflowT<T>::lc_from_rfm( rfc_counts_v &lc, rfc_value_v &level, const rfc_counts_t *rfm, rfc_flags_e flags ) const
{
    lc.resize( m_ctx.class_count );
    level.resize( m_ctx.class_count );

    return lc_from_rfm( &lc[0], &level[0], rfm, flags );
}


template< class T >
bool RainflowT<T>::lc_from_residue( rfc_counts_v &lc, rfc_value_v &level, rfc_flags_e flags ) const
{
    lc.resize( m_ctx.class_count );
    level.resize( m_ctx.class_count );

    return lc_from_residue( &lc[0], &level[0], flags );
}


template< class T >
bool RainflowT<T>::rp_get( rfc_counts_v &rp, rfc_value_v &Sa ) const
{
    rp.resize( m_ctx.class_count );
    Sa.resize( m_ctx.class_count );

    return rp_get( &rp[0], &Sa[0] );
}


template< class T >
bool RainflowT<T>::rp_from_rfm( rfc_counts_v &rp, rfc_value_v &Sa, const rfc_counts_t *rfm ) const
{
    rp.resize( m_ctx.class_count );
    Sa.resize( m_ctx.class_count );

    return rp_from_rfm( &rp[0], &Sa[0], rfm );
}


template< class T >
bool RainflowT<T>::damage_from_rp( const rfc_counts_v &counts, const rfc_value_v &Sa, double &damage, rfc_rp_damage_method_e rp_calc_type ) const
{
    return damage_from_rp( &counts[0], &Sa[0], &damage, rp_calc_type );
}


template< class T >
bool RainflowT<T>::at_init( const rfc_double_v &Sa, const rfc_double_v &Sm, 
                            double M, double Sm_rig, double R_rig, bool R_pinned, bool symmetric )
{
    if( Sa.size() != Sm.size() )
    {
        return false;
    }

    return at_init( &Sa[0], &Sm[0], (unsigned)Sa.size(), M, Sm_rig, R_rig, R_pinned, symmetric );
}


template< class T >
bool RainflowT<T>::at_init( double M, double Sm_rig, double R_rig, bool R_pinned )
{
    return at_init( /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, M, Sm_rig, R_rig, R_pinned, /*symmetric*/ false );
}


template< class T >
bool RainflowT<T>::at_transform( double Sa, double Sm, double &Sa_transformed ) const
{
    return RF::RFC_at_transform( &m_ctx, Sa, Sm, &Sa_transformed );
}


template< class T >
bool RainflowT<T>::wl_param_get( rfc_wl_param_s &wl_param ) const
{
    return RF::RFC_wl_param_get( &m_ctx, &wl_param );
}


/* Delegates */
template< class T >
bool RainflowT<T>::tp_set( size_t tp_pos, rfc_value_tuple_s *tp )
{
    if( tp_pos )
    {
        /* Alter or move existing turning point */
        if( tp_pos > m_ctx.tp_cnt )
        {
            /* Writing behind tp_cnt is not ok */
            return false;
        }

#if RFC_DH_SUPPORT
        if( tp->damage < 0.0 )
        {
            /* Negative tp->damage instructs to maintain turning points damage */
            tp->damage = m_tp[ tp_pos - 1 ].damage;
        }
#endif /*RFC_DH_SUPPORT*/

        tp->tp_pos         =  0;        /* No position information for turning points in its storage */
        m_tp[ tp_pos - 1 ] = *tp;       /* Move or replace turning point */
        tp->tp_pos         =  tp_pos;   /* Ping back the position (commonly tp lies in residue buffer) */

#if RFC_DEBUG_FLAGS
        if( m_ctx.internal.debug_flags & RF::RFC_FLAGS_LOG_WRITE_TP )
        {
            RFC_debug_fprintf( &m_ctx, stdout, 
                               "Alter tp #%lu (%g[%lu] @ %lu)\n", 
                               tp_pos, tp->value, tp->cls, tp->pos );
        }
#endif /*RFC_DEBUG_FLAGS*/

        return true;
    }
    else
    {
        /* Append (tp_pos == 0) */
        if( tp->tp_pos )
        {
#if RFC_DEBUG_FLAGS
            if( m_ctx.internal.debug_flags & RF::RFC_FLAGS_LOG_WRITE_TP )
            {
                RFC_debug_fprintf( &m_ctx, stdout, 
                                   "Alter (unchanged) tp #%lu (%g[%lu] @ %lu)\n", 
                                   tp->pos, tp->value, tp->cls, tp->pos );
            }
#endif /*RFC_DEBUG_FLAGS*/

            /* Already an element of tp stack */
            return tp->tp_pos <= m_ctx.tp_cnt;
        }
        else
        {
            size_t tp_pos = ++m_ctx.tp_cnt;

            /* Append tp at the tail */
            tp->tp_pos = 0;
            
            if( tp_pos > m_tp.size() )
            {
                m_tp.push_back( *tp );
            }
            else
            {
                m_tp[tp_pos-1] = *tp;
            }

            tp->tp_pos = tp_pos;
            m_ctx.tp_cap = m_tp.capacity();

#if RFC_DEBUG_FLAGS
            if( m_ctx.internal.debug_flags & RF::RFC_FLAGS_LOG_WRITE_TP )
            {
                RF::RFC_debug_fprintf( &m_ctx, stdout, 
                                       "Append tp #%lu (%g[%lu] @ %lu)\n", 
                                       tp->tp_pos, tp->value, tp->cls, tp->pos );
            }
#endif /*RFC_DEBUG_FLAGS*/
        }
    }

    if( m_ctx.internal.flags & RFC_FLAGS_TPAUTOPRUNE && m_ctx.tp_cnt > m_ctx.tp_prune_threshold )
    {
        return RF::RFC_tp_prune( &m_ctx, m_ctx.tp_prune_size, RF::RFC_FLAGS_TPPRUNE_PRESERVE_POS );
    }

    return true;
}


template< class T >
bool RainflowT<T>::tp_get( size_t tp_pos, rfc_value_tuple_s **tp )
{
    /* Reading behind tp_cnt is ok */
    if( !tp || !tp_pos || tp_pos > m_ctx.tp_cap )
    {
        return false;
    }

#if RFC_DEBUG_FLAGS
    if( m_ctx.internal.debug_flags & RF::RFC_FLAGS_LOG_READ_TP )
    {
        RF::RFC_debug_fprintf( &m_ctx, stdout, 
                               "Read tp #%lu (%g[%lu] @ %lu)\n", 
                               tp_pos, m_tp[tp_pos-1].value, m_tp[tp_pos-1].cls, m_tp[tp_pos-1].pos );
    }
#endif /*RFC_DEBUG_FLAGS*/

    *tp = &m_tp[tp_pos-1];

    return true;
}


template< class T >
bool RainflowT<T>::tp_inc_damage( size_t pos, double damage )
{
    if( !pos || pos > m_tp.size() )
    {
        return false;
    }

    m_tp[pos-1].damage += damage;

    return true;
}


template< class T >
void* RainflowT<T>::mem_alloc( void *ptr, size_t num, size_t size, rfc_mem_aim_e aim )
{
    if( !num || !size )
    {
        if( ptr )
        {
            FREE( ptr );
        }
        return NULL;
    }
    else
    {
        return ptr ? REALLOC( ptr, num * size ) : CALLOC( num, size );
    }
}


#ifdef RFC_TP_STORAGE

/* Define a Rainflow class with delegates for external turning point storage.
   Templates and namespaces use name mangling, which is not supported 
   for extern "C" linkage. */


/* Module static C delegates */
extern "C"
{
    static
    bool rfc_storage_tp_set( RF::rfc_ctx_s* ctx, size_t tp_pos, RF::rfc_value_tuple_s *tp )
    {
        return ctx && 
               ctx->internal.obj && 
               static_cast<Rainflow*>(ctx->internal.obj)->tp_set( tp_pos, tp );
    }

    static
    bool rfc_storage_tp_get( RF::rfc_ctx_s* ctx, size_t tp_pos, RF::rfc_value_tuple_s **tp )
    {
        return ctx && 
               ctx->internal.obj && 
               static_cast<Rainflow*>(ctx->internal.obj)->tp_get( tp_pos, tp );
    }

    static 
    bool rfc_storage_tp_inc_damage( RF::rfc_ctx_s *ctx, size_t tp_pos, double damage )
    {
        return ctx && 
               ctx->internal.obj && 
               static_cast<Rainflow*>(ctx->internal.obj)->tp_inc_damage( tp_pos, damage );
    }
}
#endif /*RFC_TP_STORAGE*/

#pragma pack(pop)
