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

namespace RF = RFC_CPP_NAMESPACE;


extern "C"
{
    static bool  tp_set           ( RF::rfc_ctx_s* ctx, size_t tp_pos, RF::rfc_value_tuple_s *tp );
    static bool  tp_get           ( RF::rfc_ctx_s* ctx, size_t tp_pos, RF::rfc_value_tuple_s **tp );
    static bool  tp_inc_damage    ( RF::rfc_ctx_s *ctx, size_t tp_pos, double damage );
    static void* mem_alloc        ( void *ptr, size_t num, size_t size, int aim );
}



class Rainflow
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
        RFC_ERROR_MEMORY                        = RF::RFC_ERROR_MEMORY,                         /**< Error on memory allocation */
        RFC_ERROR_TP                            = RF::RFC_ERROR_TP,                             /**< Error while amplitude transformation */
        RFC_ERROR_AT                            = RF::RFC_ERROR_AT,                             /**< Error while amplitude transformation */
        RFC_ERROR_LUT                           = RF::RFC_ERROR_LUT,                            /**< Error while accessing look up tables */
    };


    enum rfc_counting_method
    {
        RFC_COUNTING_METHOD_DELEGATED           = RF::RFC_COUNTING_METHOD_DELEGATED,            /**< Method must be implemented via delegator, see member cycle_find_fcn */
        RFC_COUNTING_METHOD_NONE                = RF::RFC_COUNTING_METHOD_NONE,                 /**< No counting */
        RFC_COUNTING_METHOD_4PTM                = RF::RFC_COUNTING_METHOD_4PTM,                 /**< 4 point algorithm (default) */
        RFC_COUNTING_METHOD_HCM                 = RF::RFC_COUNTING_METHOD_HCM,                  /**< 3 point algorithm, Clormann/Seeger (HCM) method */
        RFC_COUNTING_METHOD_COUNT               = RF::RFC_COUNTING_METHOD_COUNT,                /**< Number of options */
    };


    enum rfc_res_method
    {
        /* Don't change order! */
        RFC_RES_NONE                            = RF::RFC_RES_NONE,                             /**< No residual method */
        RFC_RES_IGNORE                          = RF::RFC_RES_IGNORE,                           /**< Ignore residue (same as RFC_RES_NONE) */
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

    /* Typedefs */
    typedef                 RFC_VALUE_TYPE          rfc_value_t;                                /** Input data value type */
    typedef                 RFC_COUNTS_VALUE_TYPE   rfc_counts_t;                               /** Type of counting values */
    typedef     struct      rfc_value_tuple         rfc_value_tuple_s;                          /** Tuple of value and index position */
    typedef     struct      rfc_ctx                 rfc_ctx_s;                                  /** Forward declaration (rainflow context) */
    typedef     enum        rfc_mem_aim             rfc_mem_aim_e;                              /** Memory accessing mode */
    typedef     enum        rfc_flags               rfc_flags_e;                                /** Flags, see RFC_FLAGS... */
    typedef     enum        rfc_state               rfc_state_e;                                /** Counting state, see RFC_STATE... */
    typedef     enum        rfc_error               rfc_error_e;                                /** Recent error, see RFC_ERROR... */
    typedef     enum        rfc_res_method          rfc_res_method_e;                           /** Method when count residue into matrix, see RFC_RES... */
    typedef     enum        rfc_counting_method     rfc_counting_method_e;                      /** Counting method, see RFC_COUNTING... */
    typedef     enum        rfc_rp_damage_method    rfc_rp_damage_method_e;                     /** Method when calculating damage from range pair counting, see RFC_RP_DAMAGE_CALC_METHOD... */
    typedef     enum        rfc_sd_method           rfc_sd_method_e;                            /** Spread damage method, see RFC_SD... */
    typedef     struct      rfc_class_param         rfc_class_param_s;                          /** Class parameters (width, offset, count) */
    typedef     struct      rfc_wl_param            rfc_wl_param_s;                             /** Woehler curve parameters (sd, nd, k, k2, omission) */
    typedef     struct      rfc_rfm_item            rfc_rfm_item_s;                             /** Rainflow matrix element */


    /* Memory allocation functions typedef */
    typedef     void *   ( *rfc_mem_alloc_fcn_t )   ( void *, size_t num, size_t size, rfc_mem_aim_e aim );     /** Memory allocation functor */

    /* Core functions */
    bool    init                    ( unsigned class_count, rfc_value_t class_width, rfc_value_t class_offset, 
                                      rfc_value_t hysteresis, rfc_flags_e flags );
    bool    wl_init_elementary      ( double sx, double nx, double k );
    bool    wl_init_original        ( double sd, double nd, double k );
    bool    wl_init_modified        ( double sx, double nx, double k, double k2 );
    bool    wl_init_any             ( const rfc_wl_param_s* );
    bool    clear_counts            ();
    bool    deinit                  ();
    bool    feed                    ( const rfc_value_t* data, size_t count );
    bool    cycle_process_counts    ( rfc_value_t from_val, rfc_value_t to_val, rfc_flags_e flags );
    bool    feed_scaled             ( const rfc_value_t* data, size_t count, double factor );
    bool    feed_tuple              ( rfc_value_tuple_s *data, size_t count );
    bool    finalize                ( rfc_res_method_e residual_method );
    /* Functions on rainflow matrix */           
    bool    rfm_make_symmetric      ();
    bool    rfm_get                 ( rfc_rfm_item_s **buffer, unsigned *count );
    bool    rfm_set                 ( const rfc_rfm_item_s *buffer, unsigned count, bool add_only );
    bool    rfm_peek                ( rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t *count );
    bool    rfm_poke                ( rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t count, bool add_only );
    bool    rfm_sum                 ( unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, rfc_counts_t *count );
    bool    rfm_damage              ( unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, double *damage );
    bool    rfm_check               ();
    bool    lc_get                  ( rfc_counts_t *lc, rfc_value_t *level );
    bool    lc_from_rfm             ( rfc_counts_t *lc, rfc_value_t *level, const rfc_counts_t *rfm, rfc_flags_e flags );
    bool    lc_from_residue         ( rfc_counts_t *lc, rfc_value_t *level, rfc_flags_e flags );
    bool    rp_get                  ( rfc_counts_t *rp, rfc_value_t *class_means );
    bool    rp_from_rfm             ( rfc_counts_t *rp, rfc_value_t *class_means, const rfc_counts_t *rfm );
    bool    damage_from_rp          ( const rfc_counts_t *counts, const rfc_value_t *Sa, double *damage, rfc_rp_damage_method_e rp_calc_type );
    bool    damage_from_rfm         ( const rfc_counts_t *rfm, double *damage );
    bool    wl_calc_sx              ( double s0, double n0, double k, double *sx, double nx, double  k2, double  sd, double nd );
    bool    wl_calc_sd              ( double s0, double n0, double k, double  sx, double nx, double  k2, double *sd, double nd );
    bool    wl_calc_k2              ( double s0, double n0, double k, double  sx, double nx, double *k2, double  sd, double nd );
    bool    wl_calc_sa              ( double s0, double n0, double k, double  n,  double *sa );
    bool    wl_calc_n               ( double s0, double n0, double k, double  sa, double *n );
    bool    tp_init                 ( rfc_value_tuple_s *tp, size_t tp_cap, bool is_static );
    bool    tp_init_autoprune       ( bool autoprune, size_t size, size_t threshold );
    bool    tp_prune                ( size_t count, rfc_flags_e flags );
    bool    dh_init                 ( int method, double *dh, size_t dh_cap, bool is_static );
    bool    at_init                 ( const double *Sa, const double *Sm, unsigned count, 
                                      double M, double Sm_rig, double R_rig, bool R_pinned, bool symmetric );
    bool    at_transform            ( double Sa, double Sm, double *Sa_transformed );

    /* C++ specific extensions */
    bool    feed                    ( const std::vector<rfc_value_t> data );
private:
    RF::rfc_ctx_s m_ctx;
};


bool Rainflow::init( unsigned class_count, rfc_value_t class_width, rfc_value_t class_offset, 
                     rfc_value_t hysteresis, rfc_flags_e flags = RFC_FLAGS_DEFAULT )
{
    RF::rfc_ctx_s ctx = { sizeof( RF::rfc_ctx_s ) };

    m_ctx = ctx;

    return RF::RFC_init( &m_ctx, class_count, class_width, class_offset, hysteresis, (RF::rfc_flags_e)flags );
}


bool Rainflow::wl_init_elementary( double sx, double nx, double k )
{
    return RF::RFC_wl_init_elementary( &m_ctx, sx, nx, k );
}


bool Rainflow::wl_init_original( double sd, double nd, double k )
{
    return RF::RFC_wl_init_original( &m_ctx, sd, nd, k );
}


bool Rainflow::wl_init_modified( double sx, double nx, double k, double k2 )
{
    return RF::RFC_wl_init_modified( &m_ctx, sx, nx, k, k2 );
}


bool Rainflow::wl_init_any( const rfc_wl_param_s* wl_param )
{
    return RF::RFC_wl_init_any( &m_ctx, (const RF::rfc_wl_param_s*) wl_param );
}


bool Rainflow::clear_counts()
{
    return RF::RFC_clear_counts( &m_ctx );
}


bool Rainflow::deinit()
{
    return RF::RFC_deinit( &m_ctx );
}


bool Rainflow::feed( const rfc_value_t* data, size_t count )
{
    return RF::RFC_feed( &m_ctx, (const RF::rfc_value_t*)data, count );
}


bool Rainflow::cycle_process_counts( rfc_value_t from_val, rfc_value_t to_val, rfc_flags_e flags )
{
    return RF::RFC_cycle_process_counts( &m_ctx, (RF::rfc_value_t)from_val, (RF::rfc_value_t)to_val, (RF::rfc_flags_e)flags );
}


bool Rainflow::feed_scaled( const rfc_value_t* data, size_t count, double factor )
{
    return RF::RFC_feed_scaled( &m_ctx, (const RF::rfc_value_t*)data, count, factor );
}


bool Rainflow::feed_tuple( rfc_value_tuple_s *data, size_t count )
{
    return RF::RFC_feed_tuple( &m_ctx, (RF::rfc_value_tuple_s *)data, count );
}


bool Rainflow::finalize( rfc_res_method_e residual_method = RFC_RES_IGNORE )
{
    return RF::RFC_finalize( &m_ctx, (RF::rfc_res_method_e)residual_method );
}


bool Rainflow::rfm_make_symmetric()
{
    return RF::RFC_rfm_make_symmetric( &m_ctx );
}


bool Rainflow::rfm_get( rfc_rfm_item_s **buffer, unsigned *count )
{
    return RF::RFC_rfm_get( &m_ctx, (RF::rfc_rfm_item_s **)buffer, count );
}


bool Rainflow::rfm_set( const rfc_rfm_item_s *buffer, unsigned count, bool add_only )
{
    return RF::RFC_rfm_set( &m_ctx, (const RF::rfc_rfm_item_s *)buffer, count, add_only );
}


bool Rainflow::rfm_peek( rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t *counts )
{
    return RF::RFC_rfm_peek( &m_ctx, (RF::rfc_value_t)from_val, (RF::rfc_value_t)to_val, (RF::rfc_counts_t *)counts );
}


bool Rainflow::rfm_poke( rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t counts, bool add_only )
{
    return RF::RFC_rfm_poke( &m_ctx, (RF::rfc_value_t)from_val, (RF::rfc_value_t)to_val, (RF::rfc_counts_t)counts, add_only );
}


bool Rainflow::rfm_sum( unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, rfc_counts_t *count )
{
    return RF::RFC_rfm_sum( &m_ctx, from_first, from_last, to_first, to_last, (RF::rfc_counts_t *)count );
}


bool Rainflow::rfm_damage( unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, double *damage )
{
    return RF::RFC_rfm_damage( &m_ctx, from_first, from_last, to_first, to_last, damage );
}


bool Rainflow::rfm_check()
{
    return RF::RFC_rfm_check( &m_ctx );
}


bool Rainflow::lc_get( rfc_counts_t *lc, rfc_value_t *level )
{
    return RF::RFC_lc_get( &m_ctx, (RF::rfc_counts_t *)lc, (RF::rfc_value_t *)level );
}


bool Rainflow::lc_from_rfm( rfc_counts_t *lc, rfc_value_t *level, const rfc_counts_t *rfm, rfc_flags_e flags )
{
    return RF::RFC_lc_from_rfm( &m_ctx, (RF::rfc_counts_t *)lc, (RF::rfc_value_t *)level, (const RF::rfc_counts_t *)rfm, (RF::rfc_flags_e)flags );
}


bool Rainflow::lc_from_residue( rfc_counts_t *lc, rfc_value_t *level, rfc_flags_e flags )
{
    return RF::RFC_lc_from_residue( &m_ctx, (RF::rfc_counts_t *)lc, (RF::rfc_value_t *)level, (RF::rfc_flags_e) flags );
}


bool Rainflow::rp_get( rfc_counts_t *rp, rfc_value_t *class_means )
{
    return RF::RFC_rp_get( &m_ctx, (RF::rfc_counts_t *)rp, (RF::rfc_value_t *)class_means );
}


bool Rainflow::rp_from_rfm( rfc_counts_t *rp, rfc_value_t *class_means, const rfc_counts_t *rfm )
{
    return RF::RFC_rp_from_rfm( &m_ctx, (RF::rfc_counts_t *)rp, (RF::rfc_value_t *)class_means, (const RF::rfc_counts_t *)rfm );
}


bool Rainflow::damage_from_rp( const rfc_counts_t *counts, const rfc_value_t *Sa, double *damage, rfc_rp_damage_method_e rp_calc_type )
{
    return RF::RFC_damage_from_rp( &m_ctx, (const RF::rfc_counts_t *)counts, (const RF::rfc_value_t *)Sa, damage, (RF::rfc_rp_damage_method_e)rp_calc_type );
}


bool Rainflow::damage_from_rfm( const rfc_counts_t *rfm, double *damage )
{
    return RF::RFC_damage_from_rfm( &m_ctx, (const RF::rfc_counts_t *)rfm, damage );
}


bool Rainflow::wl_calc_sx( double s0, double n0, double k, double *sx, double nx, double  k2, double  sd, double nd )
{
    RF::RFC_wl_calc_sx( &m_ctx, s0, n0, k, sx, nx, k2, sd, nd );
}


bool Rainflow::wl_calc_sd( double s0, double n0, double k, double  sx, double nx, double  k2, double *sd, double nd )
{
    return RF::RFC_wl_calc_sd( &m_ctx, s0, n0, k, sx, nx, k2, sd, nd );
}


bool Rainflow::wl_calc_k2( double s0, double n0, double k, double  sx, double nx, double *k2, double  sd, double nd )
{
    RF::RFC_wl_calc_k2( &m_ctx, s0, n0, k, sx, nx, k2, sd, nd );
}


bool Rainflow::wl_calc_sa( double s0, double n0, double k, double  n,  double *sa )
{
    return RF::RFC_wl_calc_sa( &m_ctx, s0, n0, k, n,  sa );
}


bool Rainflow::wl_calc_n( double s0, double n0, double k, double  sa, double *n )
{
    return RF::RFC_wl_calc_n( &m_ctx, s0, n0, k, sa, n );
}


bool Rainflow::tp_init( rfc_value_tuple_s *tp, size_t tp_cap, bool is_static )
{
    return RF::RFC_tp_init( &m_ctx, (RF::rfc_value_tuple_s *)tp, tp_cap, is_static );
}


bool Rainflow::tp_init_autoprune( bool autoprune, size_t size, size_t threshold )
{
    return RF::RFC_tp_init_autoprune( &m_ctx, autoprune, size, threshold );
}


bool Rainflow::tp_prune( size_t count, rfc_flags_e flags )
{
    return RF::RFC_tp_prune( &m_ctx, count, (RF::rfc_flags_e) flags );
}


bool Rainflow::dh_init( int method, double *dh, size_t dh_cap, bool is_static )
{
    return RF::RFC_dh_init( &m_ctx, method, dh, dh_cap, is_static );
}


bool Rainflow::at_init( const double *Sa, const double *Sm, unsigned count,
                        double M, double Sm_rig, double R_rig, bool R_pinned, bool symmetric )
{
    return RF::RFC_at_init( &m_ctx, Sa, Sm, count, M, Sm_rig, R_rig, R_pinned, symmetric );
}


bool Rainflow::at_transform( double Sa, double Sm, double *Sa_transformed )
{
    return RF::RFC_at_transform( &m_ctx, Sa, Sm, Sa_transformed );
}


bool Rainflow::feed( const std::vector<rfc_value_t> data )
{
    return RF::RFC_feed( &m_ctx, &data[0], data.size() );
}