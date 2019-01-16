/*
    Woehler:    (Sa/SD)^k == n/ND
    Basquin:    C         == n * Sa^b       (C = 1e22 fuer SD=1e3, ND=1e7 und b=5)

    Einfache Berechnung der BKZ (elementar):
    Eckschwingspielzahl:            ND (=1e7)
    Amplitude bei ND:               SD (=1e3)
    Neigung:                        k  (=-5)
    Signalamplitude in Klasse i:    Sa_i = ABS(Startklasse_i – Zielklasse_i) * Klassenbreite/2
    Zählung in Klasse i:            h_i
    Teilschädigung in Klasse i:     D_i = h_i/ND * (Sa_i/SD) ^ ABS(k)
                                    D_i = h_i * Sa ^ b / C
    Schadenssumme (fikt.):          BKZ = Summe( D_i )
*/



#pragma once

#undef ASSERT
#include <assert.h>
#define ASSERT(x) assert(x)

#include <vector>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include "rainflow.h"

namespace RF = RFC_CPP_NAMESPACE;

#ifndef DBL_NAN
#if _WIN32
#if _MSC_VER
#define DBL_NAN (_Nan._Double)
#endif /*_MSC_VER*/
#else /*!_WIN32*/
#define DBL_NAN NAN
#endif /*_WIN32*/
#endif /*DBL_NAN*/

#ifndef DBL_ISFINITE
#if _WIN32
#if _MSC_VER
#define DBL_ISFINITE(x) (_finite(x))
#endif /*_MSC_VER*/
#else /*!_WIN32*/
#define DBL_ISFINITE(x) (std::isfinite(x))
#endif /*_WIN32*/
#endif /*DBL_NAN*/


/** Returns the less value of \a a or \a b */
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
/** Returns the greater value of \a a or \a b */
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
/** Returns 1, if \a a >= 0, otherwise -1 */
#define SIGN2(a) ( ( (a) >= 0 ) ? 1 : -1 )
/** Returns \f$ {\mathop \rm \sgn}\left( a \right) \f$ */
#define SIGN3(a) ( ( (a) > 0 ) ? 1 : ( ( (a) < 0 ) ? -1 : 0 ) )
/** Returns \a a rounded */
#define ROUND(a) ( SIGN2(a) * floor( fabs( (double) a ) + 0.5 ) )
/** Returns \a a rounded to nearest integer */
#define IROUND(a) ( (int)SIGN2(a) * (int)( fabs( (double) a ) + 0.5 ) )
/** Returns \a a rounded to nearest positive integer */
#define NROUND(a) ( (size_t)( fabs( (double) a ) + 0.5 ) )


template<class T> class CRainflowT;
typedef CRainflowT<RF::RFC_value_type> CRainflow;

extern "C"
{
    static bool tp_set           ( RF::rfc_ctx_s* ctx, size_t tp_pos, RF::rfc_value_tuple_s *tp );
    static bool tp_get           ( RF::rfc_ctx_s* ctx, size_t tp_pos, RF::rfc_value_tuple_s **tp );
    static bool tp_inc_damage    ( RF::rfc_ctx_s *ctx, size_t tp_pos, double damage );
    static bool tp_prune         ( RF::rfc_ctx_s *ctx, size_t count, int flags );
}


template<class T>
class CRainflowT 
{
public:
    typedef enum 
    {
        RESIDUUM_NONE,             /* Kein Residuum, keine Matrix mit Residuum */
        RESIDUUM_IGNORE,           /* Residuum verwerfen */
        RESIDUUM_HALFCYCLES,       /* ASTM */
        RESIDUUM_FULLCYCLES,       /* Residuum als volle Schwingspiele werten */
        RESIDUUM_CLORMANN_SEEGER,  /* Bewertung nach Clormann-Seeger */
        RESIDUUM_REPEATED,         /* Wiederholter Durchlauf */
        RESIDUUM_RP_DIN,           /* Spannenpaar nach DIN */
    } e_residuum;

    typedef enum
    {
        SPRDAM_NONE              = -1,  /* Keine Schaedigung aufteilen */   
        SPRDAM_HALF_23           =  0,  /* Schaedigung jeweils zur Haelfte auf P2 und P3 */   
        SPRDAM_RAMP_AMPLITUDE_23 =  1,  /* Lineare Amplitude   ueber P2 bis P3 */   
        SPRDAM_RAMP_DAMAGE_23    =  2,  /* Lineare Schaedigung ueber P2 bis P3 */   
        SPRDAM_RAMP_AMPLITUDE_24 =  3,  /* Lineare Amplitude   ueber P2 bis P4 */   
        SPRDAM_RAMP_DAMAGE_24    =  4,  /* Lineare Schaedigung ueber P2 bis P4 */
        SPRDAM_FULL_P2           =  5,  /* Schaedigung auf P2 */
        SPRDAM_FULL_P3           =  6,  /* Schaedigung auf P3 */
    } e_spread_damage;

    typedef struct 
    {
        double dRange;              /* Doppelamplitude */
        double dAvrg;               /* Mittellast */
        double dCounts;             /* Schwingspiele */
    } t_stufe;

    typedef struct 
    {
        double dAmpl;               /* Amplitude */
        double dCounts;             /* Schwingspiele */
    } t_stufe_ampl;

    typedef struct 
    {
        double dLevel;              /* Klassengrenze */
        double dCounts;             /* Schwingspiele */
    } t_stufe_KGUZ;

    typedef struct 
    {
        double dFrom;               /* Startklasse */
        double dTo;                 /* Zielklasse */
        double dCounts;             /* Schwingspiele */
    } t_stufe_from_to;

    typedef struct 
    {
        int    iCno;                /* Klassennummer, base 0 */
        T      Value;               /* Messwert */
        T      ClassMean;           /* Klassenmitte */
        size_t nIdx;                /* Messwertposition, base 1 */
        size_t nIdx_TP;             /* Position in m_TurningPoints, base 1 */
    } t_res;


    typedef struct 
    {
      double ND, SD, k1, k2;
    } t_sn_curve;

    typedef RF::rfc_value_tuple_s         t_turning_point;
    typedef std::vector<double>           VectorData;
    typedef std::vector<t_res>            VectorResiduum;
    typedef std::vector<t_stufe>          VectorStufen;
    typedef std::vector<t_stufe_ampl>     VectorStufenAmpl;
    typedef std::vector<t_stufe_KGUZ>     VectorStufenKGUZ;
    typedef std::vector<t_stufe_from_to>  VectorStufenFromTo;
#if HAVE_NEOLIB
    typedef HUGEVECTOR(t_turning_point)   VectorTurningPoints;
#else
    typedef std::vector<t_turning_point>  VectorTurningPoints;
#endif

    enum 
    { 
        DEFAULT_HYSTERESIS  = -1,         /* Defaultwert wird zu jew. Klassenbreite! */
        DEFAULT_CLASS_COUNT =  100,       /* Default fuer EGDB Streckenvergleich */
        DEFAULT_DILATION    =  110,       /* Default fuer EGDB Streckenvergleich (in Prozent!)*/
        DEFAULT_WL_SD       =  1000,      /* Default fuer EGDB Streckenvergleich --> 1E3 */
        DEFAULT_WL_ND       =  10000000,  /* Default fuer EGDB Streckenvergleich --> 1E7 */
        DEFAULT_WL_k        = -5,         /* Default fuer EGDB Streckenvergleich */
        DEFAULT_ROUNDOFF    =  1000000    /* Fuer Vergleichbarkeit mit LMS TecWare */
    }; // end enum

    static const RF::RFC_counts_type HALF_CYCLE_INCREMENT = RFC_HALF_CYCLE_INCREMENT;
    static const RF::RFC_counts_type FULL_CYCLE_INCREMENT = RFC_FULL_CYCLE_INCREMENT;

protected:
    // Eingangsparameter
    T                           m_range_min, m_class_width;       /* Klassierbereichsuntergrenze, Klassenbreite */   
    int                         m_class_count;                    /* Klassenzahl */   
    T                           m_hysteresis;                     /* Rueckstellbreite */
    T                           m_dilation;                       /* Aufschlag (10% bei EGDB Streckenvergleich) */
    bool                        m_isSymmetric;                    /* Bei Symmetrie wird nur das Dreieck rechts oben der Matrix betrachtet */
    t_sn_curve                  m_SNCurve;                        /* Woehlerlinienparameter */
    size_t                      m_amount;                         /* Anzahl der zu klassierenden Samplepoints */
    int                         m_doSpreadDamage;                 /* >0 = Schaedigung ueber Teilbereich verteilen, statt nur auf die ausseren TP */

    // Ausgangswerte
    VectorTurningPoints         m_TurningPoints;                  /* Umkehrpunktfolge */
    bool                        m_isParametrized;                 /* Klassierparameter gesetzt? */
    e_residuum                  m_eResiduum;                      /* Wie wurde das Residuum bewertet */
    int                         m_iProgressState;                 /* Aktueller Fortschritt in Prozent (0-100) */
    bool                        m_bCountPending;                  /* Bei mehrfachem Aufruf von DoRainflow() true */
    size_t                      m_CurrentPos;                     /* Anzahl bisher behandelter Eingangsdaten */
    double                      m_AmplTrans_dM;                   /* Mittelspannungseinfluss */
    double                      m_AmplTrans_dR;                   /* Neues R */
    double                      m_AmplTrans_dAvrg;                /* Neuer Mittelwert */
    int                         m_AmplTrans_mode;                 /* 0 = Aus, 1 = Ziel ist R, 2 = Ziel ist Avrg */

    // Rainflow
    RF::rfc_ctx_s               m_rfc_ctx;                        /* Rainflow context */

                                CRainflowT                    ( const CRainflowT &RainflowMatrix );  /* Disable Standard Constructor */
                                CRainflowT& operator=         ( const CRainflowT &RainflowMatrix );  /* Disable Copy Constructor */

    inline T                    Min                                 ( T A, T B ) const { return ( A < B ) ? A : B; }
    inline T                    Max                                 ( T A, T B ) const { return ( A > B ) ? A : B; }
    inline T                    Sign                                ( T A ) const { return ( A < 0 ) ? -1 : 1; }
    inline T                    Abs                                 ( T A ) const { return fabs( (double) A ); }
    inline T                    Ceil                                ( T A ) const { return (int) ceil( A ); }

public:
    explicit 
    CRainflowT()
    : m_iProgressState( -1 )
    { 
        RF::rfc_ctx_s dummy = { sizeof( RF::rfc_ctx_s ) };

        m_SNCurve.SD = DEFAULT_WL_SD;
        m_SNCurve.ND = DEFAULT_WL_ND;
        m_SNCurve.k1 = DEFAULT_WL_k;
        m_SNCurve.k2 = DEFAULT_WL_k;

        m_doSpreadDamage = SPRDAM_HALF_23;  

        m_rfc_ctx                   = dummy;
        m_rfc_ctx.internal.obj      = this;
        m_rfc_ctx.tp_set_fcn        = tp_set;
        m_rfc_ctx.tp_get_fcn        = tp_get;
        m_rfc_ctx.tp_inc_damage_fcn = tp_inc_damage;
        m_rfc_ctx.tp_prune_fcn      = tp_prune;


        SetAmplTrans( DBL_NAN, DBL_NAN, false );
        ZeroInit();
    } // end of CRainflowT
    

    virtual ~CRainflowT()                                            
    { 
        RF::RFC_deinit( &m_rfc_ctx );
        ClearTurningPoints();
    } // end of ~CRainflowT
    

    // Berechnet die Klassennummer k (k>=0) fuer einen Messwert
    // bin(k) definiert die Klassengrenzen
    // bin(0) = m_range_min
    // bin(1) = m_range_min + m_class_width
    // k(x) => bin(k) <= x < bin(k+1)
    inline 
    int ClassNo( T value ) const
    {
        int cno = (int) floor( ( value - m_range_min ) / m_class_width );
        
        cno = MIN( cno, m_class_count - 1 );
        cno = MAX( cno, 0 );
        
        return cno;
    } // end of ClassNo
    

    inline 
    RF::RFC_counts_type* GetMatrix() const
    {
        return m_rfc_ctx.rfm;
    } // end of GetMatrix
    

    // Gibt Anzahl Schwingspiele x FULL_CYCLE_INCREMENT zurueck!
    int GetCountIncrements( int from, int to ) const
    {
        RF::RFC_counts_type counts;
        RF::RFC_counts_type *prf = GetMatrix();

        ASSERT( m_isParametrized );

        if( prf != NULL && m_class_count > 0 ) 
        {
            counts  = prf[ m_class_count * from + to ]; // Rainflowzaehlung
        } // end if

        return counts;
    } // end of GetCountIncrements
    

    void SetCounts( int from, int to, int counts )
    {
        ASSERT( m_isParametrized );
        ASSERT( from < m_class_count && from >= 0 && to < m_class_count && to >= 0 );

        int *prf = GetMatrix();

        if( m_isSymmetric && to < from ) 
        {
            int temp = from;
            to = from;
            from = temp;
        } // end if

        ASSERT( prf && prf[ m_class_count * from + to ] == 0 );
        prf[ m_class_count * from + to ] = counts * FULL_CYCLE_INCREMENT;
    } // end of SetCounts
    

    inline
    void SetSNCurve( double SD, double ND, double k )
    {
        if( DBL_ISFINITE( SD ) && SD > 0.0 )
        {
            m_SNCurve.SD = SD;
        } // end if
        
        if( DBL_ISFINITE( ND ) && ND > 0.0 )
        {
            m_SNCurve.ND = ND;
        } // end if
        
        if( DBL_ISFINITE( k ) && fabs( k ) >= 1.0 )
        {
            m_SNCurve.k1  = -fabs( k );
            m_SNCurve.k2  = -fabs( k );
        } // end if

        RF::RFC_wl_init_elementary( &m_rfc_ctx, SD, ND, k );
    } // end of SetSNCurve
    

    inline
    void SetSNCurve( double SD, double ND, double k1, double k2 )
    {
        if( DBL_ISFINITE( SD ) && SD > 0.0 )
        {
            m_SNCurve.SD = SD;
        } // end if
        
        if( DBL_ISFINITE( ND ) && ND > 0.0 )
        {
            m_SNCurve.ND = ND;
        } // end if
        
        if( DBL_ISFINITE( k1 ) && fabs( k1 ) >= 1.0 )
        {
            m_SNCurve.k1 = -fabs( k1 );
        } // end if
        
        if( DBL_ISFINITE( k2 ) && fabs( k2 ) >= 1.0 )
        {
            m_SNCurve.k2 = -fabs( k2 );
        }
        else
        {
            m_SNCurve.k2 = -fabs(m_SNCurve.k1);
        } // end if

        RF::RFC_wl_init_modified( &m_rfc_ctx, SD, ND, k1, k2 );
    } // end of SetSNCurve
    

    inline
    void GetSNCurve( double& SD, double& ND, double& k ) const
    {
        SD = m_SNCurve.SD;
        ND = m_SNCurve.ND;
        k  = m_SNCurve.k1;
    } // end of GetSNCurve


    inline
    void GetSNCurve( double& SD, double& ND, double& k1, double& k2 ) const
    {
        SD = m_SNCurve.SD;
        ND = m_SNCurve.ND;
        k1 = m_SNCurve.k1;
        k2 = m_SNCurve.k2;
    } // end of GetSNCurve


    inline 
    void SetAmplTrans( double dM, double dZiel, bool bZielItR )
    {
        if( !DBL_ISFINITE( dM ) || !DBL_ISFINITE( dZiel ) )
        {
            m_AmplTrans_mode  = 0;
            m_AmplTrans_dM    = DBL_NAN;
            m_AmplTrans_dR    = DBL_NAN;
            m_AmplTrans_dAvrg = DBL_NAN;

            m_rfc_ctx.at.count = 0;
        }
        else
        {
            if( bZielItR )
            {
                m_AmplTrans_mode  = 1;
                m_AmplTrans_dM    = dM;
                m_AmplTrans_dR    = dZiel;
                m_AmplTrans_dAvrg = DBL_NAN;

                RF::RFC_at_init( &m_rfc_ctx, /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, dM, 
                                             /*Sm_rig*/ 0, /*R_rig*/ dZiel, /*R_pinned*/ true, /*symmetric*/ false );
            }
            else
            {
                m_AmplTrans_mode  = 2;
                m_AmplTrans_dM    = dM;
                m_AmplTrans_dR    = DBL_NAN;
                m_AmplTrans_dAvrg = dZiel;

                RF::RFC_at_init( &m_rfc_ctx, /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, dM, 
                                             /*Sm_rig*/ dZiel, /*R_rig*/ -1, /*R_pinned*/ false, /*symmetric*/ false );
            }
        }
    }
    

    inline 
    void AddCycles( int from, int to, int cycles )
    {
        ASSERT( m_isParametrized );
        ASSERT( from < m_class_count && from >= 0 && to < m_class_count && to >= 0 );

        RF::RFC_counts_type *prf = GetMatrix();

        if( m_isSymmetric && to < from ) 
        {
            int temp = from;
            to = from;
            from = temp;
        } // end if

        ASSERT( prf );
        
        prf[ m_class_count * from + to ] += cycles * FULL_CYCLE_INCREMENT;
    } // end of AddCycles
    

    inline 
    void AddHalfCycle( int from, int to )
    {
        ASSERT( m_isParametrized );
        ASSERT( from < m_class_count && from >= 0 && to < m_class_count && to >= 0 );

        RF::RFC_counts_type *prf = GetMatrix();

        if( m_isSymmetric && to < from ) 
        {
            int temp = from;
            to = from;
            from = temp;
        } // end if

        ASSERT( prf );
        prf[ m_class_count * from + to ] += HALF_CYCLE_INCREMENT;
    } // end of AddHalfCycle
    
    
    double    NoValue                      () const { return DBL_NAN; }
    void      ClearTurningPoints           ();
    void      ZeroInit                     ();
    int       EntriesCount                 () const;
    void      MakeSymmetric                ();
    void      Parametrize                  ( double range_min,       double range_max,
                                             double range_fixed_min, double range_fixed_max, double range, 
                                             double class_count,     double class_width, 
                                             double dilation,        double hysteresis );
    int       IsParametrized               () { return m_isParametrized ? 1 : 0; }
    T         GetRangeMin                  () const { return m_range_min; }
    T         GetRangeMax                  () const { return GetRangeMin () + GetRange (); }
    T         GetRange                     () const { return m_class_count * m_class_width; }
    T         GetClassWidth                () const { return m_class_width; }
    int       GetClassCount                () const { return m_class_count; }
    T         GetHysteresis                () const { return m_hysteresis; }
    T         GetDilation                  () const { return m_dilation; }
    double    GetHysteresisToClassWidth    () const { return DBL_ISFINITE ( GetRange () ) ? ( m_hysteresis / GetRange () * GetClassCount () ) : DBL_NAN; }
    void      SetHysteresis                ( double dHysteresis ) { m_hysteresis = dHysteresis; }
    void      SetAmount                    ( size_t amount ) { m_amount = amount; }
    void      SetSpreadDamage              ( int doSpread ) { m_doSpreadDamage = doSpread; }

    void      GetStufenNewMean             ( VectorStufen &StufenOrig, VectorStufenAmpl &Stufen, double dNewMean, double dM ) const;
    void      GetStufenNewR                ( VectorStufen &StufenOrig, VectorStufenAmpl &Stufen, double dNewR, double dM ) const;

    void      DoRainflow                   ( const double *stream_in, size_t in_len, bool do_finalize = true );
    void      DoRainflow                   ( const VectorData& stream );
    template<class InputIt>
    void      DoRainflow                   ( const InputIt first, const InputIt last );
    void      CountResiduum                ( e_residuum how );
    void      GetResiduum                  ( VectorResiduum &Residuum ) const;
    void      SetResiduum                  ( const VectorResiduum &Residuum );
    void      GetStufen                    ( VectorStufen &Stufen ) const;
    void      GetStufen                    ( VectorStufenFromTo &Stufen ) const;
    void      GetStufen                    ( VectorStufenAmpl &Stufen ) const;
    void      GetStufen                    ( VectorStufenKGUZ &Stufen ) const;
    void      SetStufen                    ( const VectorStufen &Stufen );
    void      SetStufen                    ( const VectorStufenFromTo &Stufen );
    void      GetResiduumStufen            ( VectorStufen &Stufen, e_residuum how ) const;
    void      GetResiduumStufen            ( VectorStufenFromTo &Stufen, e_residuum how ) const;
    void      GetTurningPoints             ( VectorTurningPoints& TurningPoints, size_t left = 0, size_t right = 0 ) const;
    const VectorTurningPoints& 
              GetTurningPoints             () const;
    void      DetachTurningPoints          ( VectorTurningPoints &TurningPoints );
    double    CalcDamage                   ( double dAmplitude, double dAvrg, double dCounts ) const;
    double    GetBKZ                       () const;
    double    GetBKZ                       ( const VectorStufen &Stufen ) const;
    double    GetBKZ                       ( const VectorStufenAmpl &Stufen ) const;
    double    GetShapeValue                () const;
    double    GetShapeValue                ( const VectorStufenAmpl &Stufen, double &dAele ) const;
    double    GetShapeValue                ( const VectorStufenAmpl &Stufen ) const;
    double    GetSumH                      () const;
    double    GetSumH                      ( const VectorStufenAmpl &Stufen ) const;
    bool      TpSet                        ( size_t tp_pos, RF::rfc_value_tuple_s *tp );
    bool      TpGet                        ( size_t tp_pos, RF::rfc_value_tuple_s **tp );
    bool      TpIncDamage                  ( size_t tp_pos, double damage );
    bool      TpPrune                      ( size_t counts, int flags );


    template< typename ParameterMap >
    ParameterMap GetParameterMap() const
    {
        ParameterMap PMap;

        double SD, ND, k;
        
        GetSNCurve( SD, ND, k );

        PMap[ "RFM.RangeMin" ]     = GetRangeMin();
        PMap[ "RFM.RangeMax" ]     = GetRangeMax();
        PMap[ "RFM.ClassWidth" ]   = GetClassWidth();
        PMap[ "RFM.ClassCount" ]   = GetClassCount();
        PMap[ "RFM.Hysteresis" ]   = GetHysteresis();
        PMap[ "RFM.Dilation" ]     = GetDilation() - (int)GetDilation();
        PMap[ "RFM.SNCurve.SD" ]   = SD;
        PMap[ "RFM.SNCurve.ND" ]   = ND;
        PMap[ "RFM.SNCurve.k" ]    = k;
        PMap[ "RFM.BKZ" ]          = GetBKZ();
        PMap[ "RFM.H" ]            = GetSumH();
        PMap[ "RFM.v" ]            = GetShapeValue();  // Voelligkeit

        return PMap;
    } // end of GetParameterMap
}; // end of class CRainflowT


template<class T>
void CRainflowT<T>::ClearTurningPoints()
{
    m_TurningPoints.clear();
} // end of ClearTurningPoints


template<class T>
void CRainflowT<T>::ZeroInit()
{
    ClearTurningPoints();
    m_hysteresis                =  0.0;
    m_dilation                  =  0.0;
    m_range_min = m_class_width =  0.0;
    m_class_count               =  0;
    m_isSymmetric               =  false;
    m_isParametrized            =  false;
    m_eResiduum                 =  RESIDUUM_IGNORE;
    m_bCountPending             =  false;
    m_amount                    =  0;
    m_CurrentPos                =  0;

    RF::RFC_clear_counts( &m_rfc_ctx );
} // end of ZeroInit


template<class T>
int CRainflowT<T>::EntriesCount() const
{
    int i, j;
    int entries_count = 0;

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return 0;
    } // end if

    for( i = 0; i < m_class_count; i++ ) 
    {
        for( j = 0; j < m_class_count; j++ ) 
        {
            if( GetCountIncrements( i, j ) > 0 ) 
            {
                entries_count++;
            } // end if
        } /* end for */
    } /* end for */

    return entries_count;
} // end of EntriesCount


// Nur fuer Matrix nach Residuenbehandlung!
// MakeSymmetric() beeinflusst nicht die Ergebnisse aus den Zaehlverfahren:
// KGUZ und BPZ (GetStufen) zaehlen stehende und haengende Hystereseaeste!
template<class T>
void CRainflowT<T>::MakeSymmetric() 
{
    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
    } // end if

    RF::RFC_rfm_make_symmetric( &m_rfc_ctx );
} /* end of MakeSymmetric() */


template<class T>
void CRainflowT<T>::GetStufenNewR( VectorStufen &StufenOrig, VectorStufenAmpl &Stufen, 
                                   double dNewR, double dM ) const
{
    t_stufe_ampl stufe;

    Stufen.clear();

    if( RF::RFC_at_init( &m_rfc_ctx, /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, dM, 
                                     /*Sm_rig*/ 0.0, dNewR, /*R_pinned*/ true, /*symmetric*/ false ) )
    {
        for( int i = 0; i < (int)StufenOrig.size(); i++ ) 
        {
            double dTransformedAmpl;

            if( RF::RFC_at_transform( &m_rfc_ctx, StufenOrig[i].dRange / 2.0, StufenOrig[i].dAvrg, &dTransformedAmpl ) )
            {
                stufe.dAmpl   = dTransformedAmpl;
                stufe.dCounts = StufenOrig[i].dCounts;
            } else {
                stufe.dAmpl   = 0;
                stufe.dCounts = -1;  // TBD: Gibts eine bessere Loesung?
            } // end if

            Stufen.push_back( stufe );
        } // end for
    }

    RF::RFC_at_init( &m_rfc_ctx, /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, m_AmplTrans_dM, 
                                 /*Sm_rig*/ m_AmplTrans_dAvrg, /*R_rig*/ m_AmplTrans_dR, /*R_pinned*/ m_AmplTrans_mode = 2, /*symmetric*/ false );
} // end of GetStufenNewR


template<class T>
void CRainflowT<T>::GetStufenNewMean( VectorStufen &StufenOrig, VectorStufenAmpl &Stufen, 
                                      double dNewMean, double dM ) const
{
    t_stufe_ampl stufe;

    Stufen.clear();

    if( RF::RFC_at_init( &m_rfc_ctx, /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, dM, 
                                     /*Sm_rig*/ dNewMean, /*R_rig*/ -1, /*R_pinned*/ false, /*symmetric*/ false ) )
    {
        for( int i = 0; i < (int)StufenOrig.size(); i++ ) 
        {
            double dTransformedAmpl;

            if( RF::RFC_at_transform( &m_rfc_ctx, StufenOrig[i].dRange / 2.0, StufenOrig[i].dAvrg, &dTransformedAmpl ) )
            {
                stufe.dAmpl   = dTransformedAmpl;
                stufe.dCounts = StufenOrig[i].dCounts;
            } else {
                stufe.dAmpl   = 0;
                stufe.dCounts = -1;  // TBD: Gibts eine bessere Loesung?
            } // end if

            Stufen.push_back( stufe );
        } // end for
    }

    RF::RFC_at_init( &m_rfc_ctx, /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, m_AmplTrans_dM, 
                                 /*Sm_rig*/ m_AmplTrans_dAvrg, /*R_rig*/ m_AmplTrans_dR, /*R_pinned*/ m_AmplTrans_mode = 2, /*symmetric*/ false );
} // end of GetStufenNewMean


template<class T>
void CRainflowT<T>::GetStufen( VectorStufen &Stufen ) const
{
    Stufen.clear();

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
    } // end if

    if( !GetMatrix() ) 
    {
        return;
    } // end if

    for( int i = 0; i < m_class_count; i++ ) 
    {
        for( int j = 0; j < m_class_count; j++ ) 
        {
            // Haeufigkeiten durch FULL_CYCLE_INCREMENT dividieren, da DoRainflow(), und SetValue()
            // pro Schwingspiel FULL_CYCLE_INCREMENT Zaehlungen vornehmen.
            // CountResiduum() nimmt eine Zaehlung (HALF_CYCLE_INCREMENT) fuer "RESIDUUM_HALFCYCLES" vor.
            // Die Schleifen decken die volle Matrix ab!
            double counts = (double) GetCountIncrements( i, j ) / FULL_CYCLE_INCREMENT;

            if ( counts > 0.0 ) 
            {
                /* Es wird mit den Klassenmitten gerechnet */
                /* ( i + 0.5 ) - ( j + 0.5 ) --> j - i
                 * ( i + 0.5 ) + ( j + 0.5 ) --> i + j + 1
                 */
                double range = m_class_width * Abs( j - i );
                double avrg  = m_class_width * ( i + j + 1 ) / 2.0 + m_range_min;
                
                t_stufe stufe;
                stufe.dRange  = range;
                stufe.dAvrg   = avrg;
                stufe.dCounts = counts;

                Stufen.push_back( stufe );
            } // end if
        } // end for
    } // end for
} // end of GetStufen


template<class T>
void CRainflowT<T>::GetStufen( VectorStufenAmpl &Stufen ) const
{
    Stufen.clear();

    if( !m_isParametrized || !GetMatrix() ) 
    {
        ASSERT( false );
        return;
    } // end if

    for( int i = 0; i < m_class_count; i++ ) 
    {
        double counts = 0.0;

        for( int j = i; j < m_class_count; j++ ) 
        {
            // Haeufigkeiten durch FULL_CYCLE_INCREMENT dividieren, da DoRainflow(), und SetValue()
            // pro Schwingspiel FULL_CYCLE_INCREMENT Zaehlungen vornehmen.
            // CountResiduum() nimmt eine Zaehlung (HALF_CYCLE_INCREMENT) fuer "RESIDUUM_HALFCYCLES" vor.
            // Die Schleifen decken die halbe Matrix ab (Dreieck)!
            counts += (double) GetCountIncrements( j - i, j ) / FULL_CYCLE_INCREMENT;
            counts += (double) GetCountIncrements( j, j - i ) / FULL_CYCLE_INCREMENT;
        } // end for

        if( counts > 0.0 ) 
        {
            /* Es wird mit den Klassenmitten gerechnet */
            double range = m_class_width * i;
            
            t_stufe_ampl stufe;
            stufe.dAmpl   = range / 2.0;
            stufe.dCounts = counts;

            Stufen.push_back( stufe );
        } // end if
    } // end for
} // end of GetStufen


template<class T>
void CRainflowT<T>::GetStufen( VectorStufenFromTo &Stufen ) const
{
    Stufen.clear();

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
    } // end if

    if( !GetMatrix() ) 
    {
        return;
    } // end if

    for( int i = 0; i < m_class_count; i++ ) 
    {
        for( int j = 0; j < m_class_count; j++ ) 
        {
            // Schleifen decken die volle Matrix ab!
            double counts = (double) GetCountIncrements( i, j ) / FULL_CYCLE_INCREMENT;

            if ( counts > 0.0 ) 
            {
                double from = m_class_width * ( 0.5 + i ) + m_range_min;
                double   to = m_class_width * ( 0.5 + j ) + m_range_min;
                
                t_stufe_from_to stufe;
                stufe.dFrom   = from;
                stufe.dTo     = to;
                stufe.dCounts = counts;

                Stufen.push_back( stufe );
            } // end if
        } // end for
    } // end for
} // end of GetStufen


template<class T>
void CRainflowT<T>::GetStufen( VectorStufenKGUZ &Stufen ) const
{
    Stufen.clear();

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
    } // end if

    VectorStufenKGUZ KGUZ_Temp( m_class_count - 1 );
    RF::RFC_counts_type *lc     = new RF::RFC_counts_type[m_class_count];
    RF::RFC_value_type  *level  = new RF::RFC_value_type[m_class_count];

    if( RF::RFC_lc_get( &m_rfc_ctx, lc, level ) )
    {
        for( int i = 0; i < m_class_count - 1; i++ ) 
        {
            // Erster Index [0] zaehlt das Durchschreiten der ersten oberen Klassengrenze.
            KGUZ_Temp[i].dLevel  = level[i];
            KGUZ_Temp[i].dCounts = lc[i];
        } // end for
    }

    Stufen = KGUZ_Temp;

    delete[] lc;
    delete[] level;
} // end of GetStufen


template<class T>
void CRainflowT<T>::SetStufen( const VectorStufen &Stufen )
{
    size_t i;

    if( !m_isParametrized || !m_rfc_ctx.rp || !m_rfc_ctx.lc ) 
    {
        ASSERT( false );
        return;
    } // end if

    RF::RFC_clear_counts( &m_rfc_ctx );

    // Rangepair
    for( i = 0; i < Stufen.size(); i++ ) 
    {
        t_stufe stufe = Stufen[i];
        int from = ClassNo( (T)( stufe.dAvrg - stufe.dRange / 2.0 ) );
        int   to = ClassNo( (T)( stufe.dAvrg + stufe.dRange / 2.0 ) );

        ASSERT( from >= 0 && from < m_class_count 
                && to >=0 && to < m_class_count );

        AddCycles( from, to, (int)( stufe.dCounts + 0.5 ) );
    } // end for

    RF::RFC_rp_from_rfm( &m_rfc_ctx, m_rfc_ctx.rp, NULL, NULL );
    RF::RFC_lc_from_rfm( &m_rfc_ctx, m_rfc_ctx.lc, NULL, NULL, RF::RFC_FLAGS_COUNT_LC );
    RF::RFC_damage_from_rfm( &m_rfc_ctx, NULL, &m_rfc_ctx.damage );
} // end of SetStufen


template<class T>
void CRainflowT<T>::SetStufen( const VectorStufenFromTo &Stufen )
{
    size_t i;

    if( !m_isParametrized || !m_rfc_ctx.rp || !m_rfc_ctx.lc ) 
    {
        ASSERT( false );
        return;
    } // end if

    RF::RFC_clear_counts( &m_rfc_ctx );

    for( i = 0; i < Stufen.size(); i++ ) 
    {
        t_stufe_from_to stufe = Stufen[i];
        int from = stufe.dFrom;
        int   to = stufe.dTo;

        ASSERT( from >= 0 && from < m_class_count 
                && to >=0 && to < m_class_count );

        AddCycles( from, to, (int)( stufe.dCounts + 0.5 ) );
    } // end for

    RF::RFC_rp_from_rfm( &m_rfc_ctx, m_rfc_ctx.rp, NULL, NULL );
    RF::RFC_lc_from_rfm( &m_rfc_ctx, m_rfc_ctx.lc, NULL, NULL, RF::RFC_FLAGS_COUNT_LC );
    RF::RFC_damage_from_rfm( &m_rfc_ctx, NULL, &m_rfc_ctx.damage );
} // end of SetStufen


template<class T>
void CRainflowT<T>::GetResiduum( VectorResiduum &Residuum ) const
{
    Residuum.clear();

    for( int i = 0; i < m_rfc_ctx.residue_cnt; i++ )
    {
        t_res res;

        res.iCno      = m_rfc_ctx.residue[i].cls;                                   /* Klassennummer, base 0 */
        res.Value     = m_rfc_ctx.residue[i].value;                                 /* Messwert */
        res.ClassMean = m_rfc_ctx.class_width * res.iCno + m_rfc_ctx.class_offset;  /* Klassenmitte */
        res.nIdx      = m_rfc_ctx.residue[i].pos;                                   /* Messwertposition, base 1 */
        res.nIdx_TP   = m_rfc_ctx.residue[i].tp_pos;                                /* Position in m_TurningPoints, base 1 */

        Residuum.push_back( res );
    }
} // end of GetResiduum


/*
template<class T>
const typename CRainflowT<T>::VectorTurningPoints& CRainflowT<T>::GetTurningPoints() const
{
    return m_TurningPoints[0];
}


template<class T>
void CRainflowT<T>::GetTurningPoints( VectorTurningPoints &TurningPoints, size_t left, size_t right ) const
{
    if( !left && !right )
    {
        VectorTurningPoints aCopy( m_TurningPoints[0] );
        TurningPoints.swap( aCopy );
    }
    else
    {
        typename VectorTurningPoints::const_iterator it;
        size_t count = 0;

        TurningPoints.clear();

#if !HAVE_NEOLIB
        if( !m_TurningPoints[0].empty() )
        {
            if( !left )  left  = m_TurningPoints[0].front().nIdx;
            if( !right ) right = m_TurningPoints[0].back().nIdx;
        } // end if

        // Count number of values first...
        for( it = m_TurningPoints[0].begin(); it != m_TurningPoints[0].end(); it++ )
        {
            if( it->nIdx >= left && it->nIdx <= right )
            {
                count++;
            }
        } // end for

        // ...then allocate space 
        TurningPoints.reserve( count );
#endif


        for( it = m_TurningPoints[0].begin(); it != m_TurningPoints[0].end(); it++ )
        {
            if( it->nIdx >= left && it->nIdx <= right )
            {
                TurningPoints.push_back( *it );
            }
        } // end for
    } // end if
} // end of GetTurningPoints
*/

template<class T>
void CRainflowT<T>::Parametrize( double range_min,         double range_max,
                                 double range_fixed_min,   double range_fixed_max, double range_fixed, 
                                 double class_count,       double class_width, 
                                 double dilation,          double hysteresis ) 
{
    if( class_count > 0.0 ) 
    {
        class_count = floor( class_count + 0.5 );
    } else {
        class_count = DBL_NAN;
    } // end if

    int is_range_min_fixed      = DBL_ISFINITE( range_fixed_min );
    int is_range_max_fixed      = DBL_ISFINITE( range_fixed_max );
    int is_range_fixed          = DBL_ISFINITE( range_fixed );
    int is_class_count_fixed    = DBL_ISFINITE( class_count );
    int is_class_width_fixed    = DBL_ISFINITE( class_width );
    int is_dilation_fixed       = DBL_ISFINITE( dilation );

    double dilation_trunc;  // Ganzzahlanteil des Klassierbereichaufschlages
    double dilation_frac;   // Restanteil des Klassierbereichaufschlages
    
    double range = DBL_NAN;

    ZeroInit();

    if( // !DBL_ISFINITE( range_min ) || !DBL_ISFINITE( range_max ) ||
        is_range_fixed && is_range_min_fixed && is_range_max_fixed || 
        is_range_fixed && is_class_count_fixed && is_class_width_fixed ) 
    {
        ASSERT( false );
    } /* end if */

    /*
     *
     *   range
     *
     */
    if( !is_range_fixed ) 
    {
        range = static_cast<T>( class_width ) * class_count;
        
        if( !DBL_ISFINITE( range ) ) 
        {
            if( is_range_max_fixed && is_range_min_fixed ) {
                range = Abs( static_cast<T>( range_fixed_max ) - static_cast<T>( range_fixed_min ) );
            } else if( is_range_max_fixed && !is_range_min_fixed ) {
                range = Abs( static_cast<T>( range_fixed_max ) - static_cast<T>( range_min ) );
            } else if( !is_range_max_fixed && is_range_min_fixed ) {
                range = Abs( static_cast<T>( range_max ) - static_cast<T>( range_fixed_min ) );
            } else {
                range = Abs( static_cast<T>( range_max ) - static_cast<T>( range_min ) );
            } // end if
        } /* end if */
    } else {
        range = range_fixed;
    } /* end if */

    ASSERT( DBL_ISFINITE( range ) );
    
    /*
     *
     *   range_min & range_max
     *
     */
    if( !is_range_min_fixed ) 
    {
        if( is_range_max_fixed ) 
        {
            range_min = static_cast<T>( range_fixed_max ) - static_cast<T>( range );
        } else if( DBL_ISFINITE( range_max ) ) {
            range_min = static_cast<T>( range_max ) - static_cast<T>( range );
        } // end if
    } else {
        range_min = range_fixed_min;
    } /* end if */
    
    if( !is_range_max_fixed ) 
    {
        if( is_range_min_fixed ) 
        {
            range_max = static_cast<T>( range_fixed_min ) + static_cast<T>( range );
        } else if( DBL_ISFINITE( range_min ) ) {
            range_max = static_cast<T>( range_min ) + static_cast<T>( range );
        } // end if
    } else {
        range_max = range_fixed_max;
    } /* end if */

    ASSERT( DBL_ISFINITE ( range_min + range_max ) );
    ASSERT( range_max >= range_min );
    range = range_max - range_min;

    if( is_range_max_fixed && is_range_min_fixed ) 
    {
        is_range_fixed = true;
    } // end if
    
    /*
     *
     *   class_count
     *
     */
    if( is_class_width_fixed && !is_class_count_fixed ) 
    {
        ASSERT( class_width > 0.0 );
        class_count = (int)Ceil( static_cast<T>( range ) / static_cast<T>( class_width ) );
    } /* end if */
    
    if( !is_class_count_fixed && !is_class_width_fixed ) 
    {
        class_count = DEFAULT_CLASS_COUNT;
        is_class_count_fixed = 1;
    } /* end if */

    if( !is_dilation_fixed ) 
    {
        dilation = (double)DEFAULT_DILATION / 100.0;
        is_dilation_fixed = 1;
    } // end if
    
    dilation       = fabs( dilation );
    dilation_trunc = (int)(dilation + 1e-4);     // Ganzzahlanteil des Klassierbereichsaufschlages
    dilation_frac  = dilation - dilation_trunc;  // Restanteil des Klassierbereichsaufschlages
    dilation_frac  = ( dilation_frac < 0.0 ) ? 0.0 : dilation_frac;
    
    dilation = dilation_trunc + dilation_frac;
    
    //ASSERT( class_count > 2 );
    ASSERT( DBL_ISFINITE( dilation ) && dilation >= 0.0 );

    /* Defaultwerte aus dem LMS TecWare "Streckenvergleich"
     * - 100 Klassen
     * - Range (Max-Min) + Dilation  (--> 10%)
     * - Hysterese = 1x Klassenbreite
     */
    

    /* Algorithmus LMS TecWare:
        range_min = min(Originaldaten)
        range_max = max(Originaldaten)
        range = range_max - range_min
        class_width = range / ( class_count - 1 ) 
        range_min = range_min - class_width / 2.0 
        range_max = range_max + class_width / 2.0 
        range = range_max - range_min 

                                v--- wo kommt dilation her?
                                     DilationPercent ist Parameter (z.B. 10 %)
                                     dilation = abs(range_max - range_min)  * DilationPercent / 100.0

        range_min = range_min - dilation / 2.0  
        range_max = range_max + dilation / 2.0 

        range = range_max - range_min

        range_min = floor( range_min * 1E6 ) / 1E6 
        range_max = floor( ( range_min + range ) * 1E6 ) / 1E6      
                                               ^-- simuliert float Rechnung, passiert bei uns aber 
                                                   implizit und mit 6.5 signifikanten Stellen, das 
                                                   koennte man mit 3162277 statt 1E6 probieren...
        range = range_max - range_min 
        class_width = range / class_count

        In der Klassierung werden die Punkte quantifiziert (Klassen 1 bis 100):
        class_no = floor( ( value - class_min ) / class_with ) + 1
    */

    if( range < 1E-5 * 2 * Max( 1, Max( Abs( range_min ), Abs( range_max ) ) ) ) 
    {
        range = 1E-5 * 2 * Max( 1, Max( Abs( range_min ), Abs( range_max ) ) );
    } /* end if */
    
    
    /* Einfache Klassenbreite, 2x halbe Klassenbreite Zuschlag, wenn dilation_trunc==1 */
    if( !is_class_width_fixed ) class_width = static_cast<T>( range ) / IROUND( class_count - dilation_trunc );
    /* Halbe Klassenbreite auch an den Klassiergrenzen aufschlagen */
    if( !is_range_min_fixed )   range_min  -= static_cast<T>( class_width ) * dilation_trunc / 2.0f;
    if( !is_range_max_fixed )   range_max  += static_cast<T>( class_width ) * dilation_trunc / 2.0f;
    /* Klassierbereich neu bestimmen */
    if( !is_range_fixed)        range       = static_cast<T>( range_max ) - static_cast<T>( range_min );

    // Dilation unabhaengig von fixierten Werten immer anwenden (muss auf 0 gesetzt werden, wenn nicht gewuenscht!)
    range_min -= static_cast<T>( range )     * static_cast<T>( dilation_frac ) / 2.0f;
    range_max += static_cast<T>( range )     * static_cast<T>( dilation_frac ) / 2.0f;
    range      = static_cast<T>( range_max ) - static_cast<T>( range_min );

    
    if( DEFAULT_ROUNDOFF > 0 ) 
    {
        range_min = static_cast<T>( floor( static_cast<T>( range_min ) * DEFAULT_ROUNDOFF ) / DEFAULT_ROUNDOFF );
        range_max = static_cast<T>( floor( static_cast<T>( range_max ) * DEFAULT_ROUNDOFF ) / DEFAULT_ROUNDOFF );
        range = static_cast<T>( range_max ) - static_cast<T>( range_min );
    } // end if

    if( !is_class_width_fixed ) class_width = static_cast<T>( range ) / IROUND( class_count );
    if( !is_class_count_fixed ) class_count = (int)Ceil( static_cast<T>( range ) / static_cast<T>( class_width ) );

    if( !DBL_ISFINITE( hysteresis ) || hysteresis == DEFAULT_HYSTERESIS ) 
    {
        hysteresis = class_width;
    } /* end if */

    ASSERT( DBL_ISFINITE( hysteresis ) && hysteresis >= 0.0 );

    m_class_count = IROUND( class_count );
    m_class_width = static_cast<T>( class_width );
    m_hysteresis  = static_cast<T>( hysteresis );
    m_range_min   = static_cast<T>( range_min );
    m_dilation    = static_cast<T>( dilation );

    if( m_class_width <= 0.0 )
    {
        m_class_width = 10.0 * DBL_EPSILON;
    }

    m_isParametrized = ( m_class_count > 0 ) && ( m_class_width > 0.0 );
} // end of Parametrize


template<class T>
void CRainflowT<T>::DoRainflow( const double *stream_in, size_t in_len, bool do_finalize ) 
{
    const   double*             pin;                        /* Pointer auf Zeitreihe (Eingangsdaten) */
            double              dummy = 0.0;
            size_t              pos;                        /* Aktuelle Position in der Zeitreihe (pin) */
            T                   d0, dl;                     /* Rueckwaertsgerichtete Gradienten ueber die letzten 3 Punkte (res[end-1],res[end],*pin) */

    if( !m_isParametrized || ( !stream_in && in_len ) ) 
    {
        ASSERT( false );
        return;
    } // end if

    ASSERT( m_class_width > 0.0 && m_class_count > 1 );

    pin = stream_in ? stream_in : &dummy;

    if( !m_bCountPending ) 
    {
        m_CurrentPos = 0;

        //ReDimMatrix();             /* Neue Matrix anlegen */
        //ClearResiduum(0);          /* Altes Residuum loeschen */
        //ClearTurningPoints(0);     /* Umkehrpunkte loeschen */
        
        m_bCountPending = !do_finalize;
    }
    
    ASSERT( pin );

    /* Schleife ueber alle Punkte */
    RF::RFC_feed( &m_rfc_ctx, pin, in_len );

    if( do_finalize )
    {
        switch( m_eResiduum )
        {
            case RESIDUUM_NONE:
                RF::RFC_finalize( &m_rfc_ctx, RF::rfc_ctx_s::RFC_RES_NONE ); break;
            case RESIDUUM_IGNORE:
                RF::RFC_finalize( &m_rfc_ctx, RF::rfc_ctx_s::RFC_RES_IGNORE ); break;
            case RESIDUUM_HALFCYCLES:
                RF::RFC_finalize( &m_rfc_ctx, RF::rfc_ctx_s::RFC_RES_HALFCYCLES ); break;
            case RESIDUUM_FULLCYCLES:
                RF::RFC_finalize( &m_rfc_ctx, RF::rfc_ctx_s::RFC_RES_FULLCYCLES ); break;
            case RESIDUUM_REPEATED:
                RF::RFC_finalize( &m_rfc_ctx, RF::rfc_ctx_s::RFC_RES_REPEATED ); break;
            case RESIDUUM_CLORMANN_SEEGER:
                RF::RFC_finalize( &m_rfc_ctx, RF::rfc_ctx_s::RFC_RES_CLORMANN_SEEGER ); break;
            case RESIDUUM_RP_DIN:
                RF::RFC_finalize( &m_rfc_ctx, RF::rfc_ctx_s::RFC_RES_RP_DIN45667 ); break;
            default:
                ASSERT( false );
                return;
        }
    }
}


template<class T>
void CRainflowT<T>::DoRainflow( const VectorData& stream ) 
{
    if( !stream.empty () ) 
    {
        DoRainflow( &stream[0], stream.size() );
    }
    else
    {
        double dDummy = 0.0;

        DoRainflow( &dDummy, /*in_len*/ 0, true );
    } // end if
} // end of DoRainflow


template<class T>
template<class InputIt>
void CRainflowT<T>::DoRainflow( InputIt first, InputIt last )
{
    size_t distance = 0;
    std::vector<T> stride( 1024 );
    
    if( first == last )
    {
        double dDummy = 0.0;
        DoRainflow( &dDummy, /*in_len*/ 0, /*do_finalize*/ true );
    } else {
        for( InputIt it = first; it != last; it += distance )
        {
            distance = std::min<size_t>( std::distance( it, last ), stride.size() );
            std::copy( it, it + distance, stride.begin() );
            DoRainflow( &stride[0], distance, /*do_finalize*/ it + distance == last );
        }
    }
}


template<class T> 
double CRainflowT<T>::CalcDamage( double dAmplitude, double dAvrg, double dCounts ) const
{
    double dDamage;
    double dLW;

    RF::RFC_at_transform( (void*)&m_rfc_ctx, dAmplitude, dAvrg, &dAmplitude );

    if( dAmplitude > m_SNCurve.SD )
    {
        dLW = pow( dAmplitude / m_SNCurve.SD, m_SNCurve.k1 ) * m_SNCurve.ND;
    }
    else
    {
        dLW = pow( dAmplitude / m_SNCurve.SD, m_SNCurve.k2 ) * m_SNCurve.ND;
    }

    dDamage = dCounts / dLW;
    
    return dDamage;
} // end of CalcDamage


template<class T>
double CRainflowT<T>::GetBKZ( const VectorStufen& Stufen ) const 
{
    double dDamageSum = 0.0;

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    } // end if

    for( int i = 0; i < (int)Stufen.size(); i++ ) 
    {
        dDamageSum += CalcDamage( Stufen[i].dRange / 2.0, Stufen[i].dAvrg, Stufen[i].dCounts );
    } // end for

    return dDamageSum;
} // end of GetBKZ


template<class T>
double CRainflowT<T>::GetBKZ( const VectorStufenAmpl& Stufen ) const 
{
    double dDamageSum = 0.0;

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    } // end if

    for( int i = 0; i < (int)Stufen.size(); i++ ) 
    {
        dDamageSum += CalcDamage( Stufen[i].dAmpl, 0.0 /* dAvrg */, Stufen[i].dCounts );
    } // end for

    return dDamageSum;
} // end of GetBKZ


template<class T>
double CRainflowT<T>::GetBKZ() const 
{
    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    } // end if

    VectorStufen Stufen;
    GetStufen( Stufen );

    return GetBKZ( Stufen );
} // end of GetBKZ


template<class T>
double CRainflowT<T>::GetSumH( const VectorStufenAmpl& Stufen ) const 
{
    double dAmount = 0.0;

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    } // end if

    for( int i = 0; i < (int)Stufen.size(); i++ ) 
    {
        dAmount += Stufen[i].dCounts;
    } // end for

    return dAmount;
} // end of GetAmount


template<class T>
double CRainflowT<T>::GetSumH() const 
{
    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    } // end if

    VectorStufenAmpl Stufen;
    GetStufen( Stufen );

    return GetSumH( Stufen );
} // end of GetAmount


template<class T>
double CRainflowT<T>::GetShapeValue( const VectorStufenAmpl& Stufen, double& dAele ) const 
{
    double dV       = 0.0;
    double dMaxAmpl = 0;    // Sa,1
    double dSumH    = GetSumH( Stufen );
    double dSNk1    = fabs( m_SNCurve.k1 );
    double dSNk2    = fabs( m_SNCurve.k2 );
    double dSD      = m_SNCurve.SD;

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    } // end if

    // Maximale Amplitude im Kollektiv suchen
    for( int i = 0; i < (int)Stufen.size(); i++ ) 
    {
        if( Stufen[i].dCounts > 0 && Stufen[i].dAmpl > dMaxAmpl )
        {
            dMaxAmpl = Stufen[i].dAmpl;
        } // end if
    } // end for

    // Voelligkeitsgrad berechnen 
    // (Betriebsfestigkeit, Verfahren und Daten zur Bauteilberechnung, Erwin Haibach (2006), Springer Verlag)
    for( int i = 0; i < (int)Stufen.size(); i++ ) 
    {
        if( Stufen[i].dAmpl > dSD )
        {
            dV += Stufen[i].dCounts / dSumH * pow( Stufen[i].dAmpl / dMaxAmpl, dSNk1 );
        }
        else
        {
            // dV += Stufen[i].dCounts / dSumH * pow( Stufen[i].dAmpl / dMaxAmpl, dSNk2 );  // TODO
            dV += Stufen[i].dCounts / dSumH * pow( Stufen[i].dAmpl / dMaxAmpl, dSNk1 );
        }
    } // end for

    dAele = 1.0 / dV;  // Abstand zwischen Bauteil- und Lebensdauerlinie, Miner-elementar

    // Voelligkeitsmass zurueckgeben
    return pow( dV, 1.0 / dSNk1 );  // TODO
} // end of GetShapeValue


template<class T>
double CRainflowT<T>::GetShapeValue( const VectorStufenAmpl& Stufen ) const 
{
    double dAele;

    return GetShapeValue( Stufen, dAele );
}


template<class T>
double CRainflowT<T>::GetShapeValue() const 
{
    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    } // end if

    VectorStufenAmpl Stufen;
    GetStufen( Stufen );
    
    return GetShapeValue( Stufen );
} // end of GetShapeValue


template<class T>
bool CRainflowT<T>::TpSet( size_t tp_pos, RF::rfc_value_tuple_s *tp )
{
    m_TurningPoints.push_back( *tp );
    tp->tp_pos = m_TurningPoints.size();

    return true;
}


template<class T>
bool CRainflowT<T>::TpGet( size_t tp_pos, RF::rfc_value_tuple_s **tp )
{
    ASSERT( tp_pos <= m_TurningPoints.size() );

    *tp = &m_TurningPoints[tp_pos-1];

    return true;
}


template<class T>
bool CRainflowT<T>::TpIncDamage( size_t tp_pos, double damage )
{
    ASSERT( tp_pos <= m_TurningPoints.size() );

    m_TurningPoints[tp_pos-1].damage += damage;

    return true;
}


template<class T>
bool CRainflowT<T>::TpPrune( size_t counts, int flags )
{
    ASSERT( false );
    return true;
}


extern "C" {

static
bool tp_set( RF::rfc_ctx_s* ctx, size_t tp_pos, RF::rfc_value_tuple_s *tp )
{
    ASSERT( !tp_pos && tp );

    if( !tp->tp_pos )
    {
        CRainflow *obj;

        obj = static_cast<CRainflow*>( ctx->internal.obj );
        return obj->TpSet( tp_pos, tp );
    }

    return false;
}

static
bool tp_get( RF::rfc_ctx_s* ctx, size_t tp_pos, RF::rfc_value_tuple_s **tp )
{
    ASSERT( tp_pos && tp );

    if( tp_pos )
    {
        CRainflow *obj;

        obj = static_cast<CRainflow*>( ctx->internal.obj );
        return obj->TpGet( tp_pos, tp );
    }

    return false;
}

static
bool tp_inc_damage( RF::rfc_ctx_s *ctx, size_t tp_pos, double damage )
{
    ASSERT( tp_pos );

    if( tp_pos )
    {
        CRainflow *obj;

        obj = static_cast<CRainflow*>( ctx->internal.obj );
        return obj->TpIncDamage( tp_pos, damage );
    }

    return false;
}

static
bool tp_prune( RF::rfc_ctx_s *ctx, size_t count, int flags )
{
    CRainflow *obj;

    ASSERT( false );

    obj = static_cast<CRainflow*>( ctx->internal.obj );
    return obj->TpPrune( count, flags );
}

} // extern "C"
