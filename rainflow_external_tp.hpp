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

#if HAVE_NEOLIB
    #define RFC_TP_STORAGE      HUGEVECTOR(RF::rfc_value_tuple_s)
#else
    #define RFC_TP_STORAGE      std::vector<RF::rfc_value_tuple_s>
#endif
#include "rainflow.hpp"

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


class CRainflow : protected Rainflow
{
public:
    typedef enum 
    {
        RESIDUUM_NONE                   = RFC_RES_NONE,             /* Kein Residuum, keine Matrix mit Residuum */
        RESIDUUM_IGNORE                 = RFC_RES_IGNORE,           /* Residuum verwerfen */
        RESIDUUM_HALFCYCLES             = RFC_RES_HALFCYCLES,       /* ASTM */
        RESIDUUM_FULLCYCLES             = RFC_RES_FULLCYCLES,       /* Residuum als volle Schwingspiele werten */
        RESIDUUM_CLORMANN_SEEGER        = RFC_RES_CLORMANN_SEEGER,  /* Bewertung nach Clormann-Seeger */
        RESIDUUM_REPEATED               = RFC_RES_REPEATED,         /* Wiederholter Durchlauf */
        RESIDUUM_RP_DIN                 = RFC_RES_RP_DIN45667,      /* Spannenpaar nach DIN */
    } e_residuum;

    typedef enum
    {
        SPRDAM_NONE              = RFC_SD_NONE,                 /* Keine Schaedigung aufteilen */   
        SPRDAM_HALF_23           = RFC_SD_HALF_23,              /* Schaedigung jeweils zur Haelfte auf P2 und P3 */   
        SPRDAM_RAMP_AMPLITUDE_23 = RFC_SD_RAMP_AMPLITUDE_23,    /* Lineare Amplitude   ueber P2 bis P3 */   
        SPRDAM_RAMP_DAMAGE_23    = RFC_SD_RAMP_DAMAGE_23,       /* Lineare Schaedigung ueber P2 bis P3 */   
        SPRDAM_RAMP_AMPLITUDE_24 = RFC_SD_RAMP_AMPLITUDE_24,    /* Lineare Amplitude   ueber P2 bis P4 */   
        SPRDAM_RAMP_DAMAGE_24    = RFC_SD_RAMP_DAMAGE_24,       /* Lineare Schaedigung ueber P2 bis P4 */
        SPRDAM_FULL_P2           = RFC_SD_FULL_P2,              /* Schaedigung auf P2 */
        SPRDAM_FULL_P3           = RFC_SD_FULL_P3,              /* Schaedigung auf P3 */
    } e_spread_damage;

    typedef struct 
    {
        double dRange;              /* Doppelamplitude */
        double dAvrg;               /* Mittellast */
        double dCounts;             /* Anzahl Schwingspiele */
    } t_stufe;

    typedef struct 
    {
        double dAmpl;               /* Amplitude */
        double dCounts;             /* Anzahl Schwingspiele */
    } t_stufe_ampl;

    typedef struct 
    {
        double dLevel;              /* Klassengrenze */
        double dCounts;             /* Anzahl Schwingspiele */
    } t_stufe_KGUZ;

    typedef struct 
    {
        double dFrom;               /* Startklasse */
        double dTo;                 /* Zielklasse */
        double dCounts;             /* Anzahl Schwingspiele */
    } t_stufe_from_to;

    typedef struct 
    {
        rfc_value_t     Value;               /* Messwert */
        int             iCno;                /* Klassennummer, base 0 */
        double          ClassMean;           /* Klassenmitte */
        size_t          nIdx;                /* Messwertposition, base 1 */
        size_t          nIdx_TP;             /* Position in m_TurningPoints, base 1 */
    } t_res;


    typedef struct 
    {
      double ND, SD, k1, k2;
    } t_sn_curve;

    typedef rfc_tp_storage                VectorTurningPoints;
    typedef rfc_value_tuple_s             t_turning_point;
    typedef rfc_value_v                   VectorData;
    typedef std::vector<t_res>            VectorResiduum;
    typedef std::vector<t_stufe>          VectorStufen;
    typedef std::vector<t_stufe_ampl>     VectorStufenAmpl;
    typedef std::vector<t_stufe_KGUZ>     VectorStufenKGUZ;
    typedef std::vector<t_stufe_from_to>  VectorStufenFromTo;


    enum 
    { 
        DEFAULT_HYSTERESIS  = -1,         /* Defaultwert wird zu jew. Klassenbreite! */
        DEFAULT_CLASS_COUNT =  100,       /* Default fuer EGDB Streckenvergleich */
        DEFAULT_DILATION    =  110,       /* Default fuer EGDB Streckenvergleich (in Prozent!)*/
        DEFAULT_WL_SD       =  1000,      /* Default fuer EGDB Streckenvergleich --> 1E3 */
        DEFAULT_WL_ND       =  10000000,  /* Default fuer EGDB Streckenvergleich --> 1E7 */
        DEFAULT_WL_k        = -5,         /* Default fuer EGDB Streckenvergleich */
        DEFAULT_ROUNDOFF    =  1000000    /* Fuer Vergleichbarkeit mit LMS TecWare */
    };

protected:
    // Eingangsparameter
    rfc_value_t                 m_range_min, m_class_width;       /* Klassierbereichsuntergrenze, Klassenbreite */   
    int                         m_class_count;                    /* Klassenzahl */   
    rfc_value_t                 m_hysteresis;                     /* Rueckstellbreite */
    rfc_value_t                 m_dilation;                       /* Aufschlag (10% bei EGDB Streckenvergleich) */
    bool                        m_isSymmetric;                    /* Bei Symmetrie wird nur das Dreieck rechts oben der Matrix betrachtet */
    t_sn_curve                  m_SNCurve;                        /* Woehlerlinienparameter */
    size_t                      m_amount;                         /* Anzahl der zu klassierenden Samplepoints */
    int                         m_doSpreadDamage;                 /* >0 = Schaedigung ueber Teilbereich verteilen, statt nur auf die ausseren TP */

    // Ausgangswerte
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
                                CRainflow                           ( const CRainflow &other );  /* Disable Standard Constructor */
                                CRainflow& operator=                ( const CRainflow &other );  /* Disable Copy Constructor */

    inline rfc_value_t          Min                                 ( rfc_value_t A, rfc_value_t B ) const { return ( A < B ) ? A : B; }
    inline rfc_value_t          Max                                 ( rfc_value_t A, rfc_value_t B ) const { return ( A > B ) ? A : B; }
    inline rfc_value_t          Sign                                ( rfc_value_t A ) const { return ( A < 0 ) ? -1 : 1; }
    inline rfc_value_t          Abs                                 ( rfc_value_t A ) const { return fabs( (double) A ); }
    inline rfc_value_t          Ceil                                ( rfc_value_t A ) const { return (int) ceil( A ); }

public:
    explicit 
    CRainflow()
    : m_iProgressState( -1 ),
      m_doSpreadDamage( SPRDAM_HALF_23 )
    { 
        Rainflow::init( /*class_count*/ 0, /*class_width*/ 0, /*class_offset*/ 0, /*hysteresis*/ 0 );
        
        SetSNCurve( DEFAULT_WL_SD, DEFAULT_WL_ND, DEFAULT_WL_k );
        SetAmplTrans( DBL_NAN, DBL_NAN, false );
        ZeroInit();
    }
    

    virtual ~CRainflow()
    { 
        Rainflow::deinit();
        ClearTurningPoints();
    }
    

    // Berechnet die Klassennummer k (k>=0) fuer einen Messwert
    // bin(k) definiert die Klassengrenzen
    // bin(0) = m_range_min
    // bin(1) = m_range_min + m_class_width
    // k(x) => bin(k) <= x < bin(k+1)
    inline 
    int ClassNo( rfc_value_t value ) const
    {
        int cno = (int) floor( ( value - m_range_min ) / m_class_width );
        
        cno = MIN( cno, m_class_count - 1 );
        cno = MAX( cno, 0 );
        
        return cno;
    }


    inline
    double ClassMean( unsigned class_no ) const
    {
        return m_class_width * ( 0.5 + class_no ) + m_range_min;
    }
    

    inline 
    const rfc_counts_t* GetMatrix() const
    {
        return Rainflow::rfm_storage();
    }
    

    inline 
    rfc_counts_t* GetMatrix()
    {
        return Rainflow::rfm_storage();
    }
    

    // Gibt Anzahl Schwingspiele x .full_inc zurueck!
    // from und to sind base 0
    int GetCountIncrements( int from, int to ) const
    {
        rfc_counts_t counts;

        if( !Rainflow::rfm_peek( from, to, &counts ) )
        {
            ASSERT( false );
            counts = -1;
        }

        return (int)counts;
    }
    

    void SetCounts( int from, int to, int counts )
    {
        if( m_isSymmetric && to < from ) 
        {
            int temp = from;
            to = from;
            from = temp;
        }

        if( !Rainflow::rfm_poke( from, to, full_inc() * counts, /*add_only*/ false ) )
        {
            ASSERT( false );
        }
    }
    

    inline
    void SetSNCurve( double SD, double ND, double k )
    {
        if( DBL_ISFINITE( SD ) && SD > 0.0 )
        {
            m_SNCurve.SD = SD;
        }
        
        if( DBL_ISFINITE( ND ) && ND > 0.0 )
        {
            m_SNCurve.ND = ND;
        }
        
        if( DBL_ISFINITE( k ) && fabs( k ) >= 1.0 )
        {
            m_SNCurve.k1  = -fabs( k );
            m_SNCurve.k2  = -fabs( k );
        }

        if( !Rainflow::wl_init_elementary( SD, ND, k ) )
        {
            ASSERT( false );
        }
    }
    

    inline
    void SetSNCurve( double SD, double ND, double k1, double k2 )
    {
        if( DBL_ISFINITE( SD ) && SD > 0.0 )
        {
            m_SNCurve.SD = SD;
        }
        
        if( DBL_ISFINITE( ND ) && ND > 0.0 )
        {
            m_SNCurve.ND = ND;
        }
        
        if( DBL_ISFINITE( k1 ) && fabs( k1 ) >= 1.0 )
        {
            m_SNCurve.k1 = -fabs( k1 );
        }
        
        if( DBL_ISFINITE( k2 ) && fabs( k2 ) >= 1.0 )
        {
            m_SNCurve.k2 = -fabs( k2 );
        }
        else
        {
            m_SNCurve.k2 = -fabs(m_SNCurve.k1);
        }

        if( !Rainflow::wl_init_modified( m_SNCurve.SD, m_SNCurve.ND, m_SNCurve.k1, m_SNCurve.k2 ) )
        {
            ASSERT( false );
        }
    }
    

    inline
    void GetSNCurve( double& SD, double& ND, double& k ) const
    {
        SD = m_SNCurve.SD;
        ND = m_SNCurve.ND;
        k  = m_SNCurve.k1;
    }


    inline
    void GetSNCurve( double& SD, double& ND, double& k1, double& k2 ) const
    {
        SD = m_SNCurve.SD;
        ND = m_SNCurve.ND;
        k1 = m_SNCurve.k1;
        k2 = m_SNCurve.k2;
    }


    inline 
    void SetAmplTrans( double dM, double dZiel, bool bZielItR )
    {
        if( !DBL_ISFINITE( dM ) || !DBL_ISFINITE( dZiel ) )
        {
            m_AmplTrans_mode  = 0;
            m_AmplTrans_dM    = DBL_NAN;
            m_AmplTrans_dR    = DBL_NAN;
            m_AmplTrans_dAvrg = DBL_NAN;

            m_ctx.at.count    = 0;
        }
        else
        {
            if( bZielItR )
            {
                m_AmplTrans_mode  = 1;
                m_AmplTrans_dM    = dM;
                m_AmplTrans_dR    = dZiel;
                m_AmplTrans_dAvrg = DBL_NAN;

                if( !Rainflow::at_init( /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, dM, 
                                        /*Sm_rig*/ 0, /*R_rig*/ dZiel, /*R_pinned*/ true, /*symmetric*/ false ) )
                {
                    ASSERT( false );
                }
            }
            else
            {
                m_AmplTrans_mode  = 2;
                m_AmplTrans_dM    = dM;
                m_AmplTrans_dR    = DBL_NAN;
                m_AmplTrans_dAvrg = dZiel;

                if( !Rainflow::at_init( /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, dM, 
                                        /*Sm_rig*/ dZiel, /*R_rig*/ -1, /*R_pinned*/ false, /*symmetric*/ false ) )
                {
                    ASSERT( false );
                }
            }
        }
    }
    

    inline 
    void AddCycles( int from, int to, int cycles )
    {
        if( m_isSymmetric && to < from ) 
        {
            int temp = from;
            to = from;
            from = temp;
        }

        if( !Rainflow::rfm_poke( from, to, full_inc() * cycles, /*add_only*/ true ) )
        {
            ASSERT( false );
        }
    }
    

    inline 
    void AddHalfCycle( int from, int to )
    {
        if( m_isSymmetric && to < from ) 
        {
            int temp = from;
            to = from;
            from = temp;
        }

        if( !Rainflow::rfm_poke( from, to, m_ctx.half_inc, /*add_only*/ true ) )
        {
            ASSERT( false );
        }
    }
    
    
    double              NoValue                      () const { return DBL_NAN; }
    void                ClearTurningPoints           ();
    void                ZeroInit                     ();
    int                 EntriesCount                 () const;
    void                MakeSymmetric                ();
    void                Parametrize                  ( double range_min,       double range_max,
                                                       double range_fixed_min, double range_fixed_max, double range, 
                                                       double class_count,     double class_width, 
                                                       double dilation,        double hysteresis );
    int                 IsParametrized               () { return m_isParametrized ? 1 : 0; }
    rfc_value_t         GetRangeMin                  () const { return m_range_min; }
    rfc_value_t         GetRangeMax                  () const { return GetRangeMin () + GetRange (); }
    rfc_value_t         GetRange                     () const { return m_class_count * m_class_width; }
    rfc_value_t         GetClassWidth                () const { return m_class_width; }
    int                 GetClassCount                () const { return m_class_count; }
    rfc_value_t         GetHysteresis                () const { return m_hysteresis; }
    rfc_value_t         GetDilation                  () const { return m_dilation; }
    double              GetHysteresisToClassWidth    () const { return DBL_ISFINITE ( GetRange () ) ? ( m_hysteresis / GetRange () * GetClassCount () ) : DBL_NAN; }
    void                SetHysteresis                ( double dHysteresis ) { m_hysteresis = dHysteresis; }
    void                SetAmount                    ( size_t amount ) { m_amount = amount; }
    void                SetSpreadDamage              ( int doSpread ) { m_doSpreadDamage = doSpread; }
    void                SetResiduumType              ( e_residuum how ) { m_eResiduum = how; }

    void                GetStufenNewMean             ( VectorStufen &StufenOrig, VectorStufenAmpl &Stufen, double dNewMean, double dM );
    void                GetStufenNewR                ( VectorStufen &StufenOrig, VectorStufenAmpl &Stufen, double dNewR, double dM );

    void                DoRainflow                   ( const double *stream_in, size_t in_len, bool do_finalize = true );
    void                DoRainflow                   ( const VectorData& stream );
    template<class InputIt>
    void                DoRainflow                   ( const InputIt first, const InputIt last );
    void                Finalize                     ( e_residuum how );
    void                GetResiduum                  ( VectorResiduum &Residuum ) const;
    void                SetResiduum                  ( const VectorResiduum &Residuum );
    void                GetStufen                    ( VectorStufen &Stufen ) const;
    void                GetStufen                    ( VectorStufenFromTo &Stufen ) const;
    void                GetStufen                    ( VectorStufenAmpl &Stufen ) const;
    void                GetStufen                    ( VectorStufenKGUZ &Stufen ) const;
    void                SetStufen                    ( const VectorStufen &Stufen );
    void                SetStufen                    ( const VectorStufenFromTo &Stufen );
    void                GetResiduumStufen            ( VectorStufen &Stufen, rfc_res_method how ) const;
    void                GetResiduumStufen            ( VectorStufenFromTo &Stufen, rfc_res_method how ) const;
    void                GetTurningPoints             ( VectorTurningPoints& TurningPoints, size_t left = 0, size_t right = 0 ) const;
    const VectorTurningPoints& 
                        GetTurningPoints             () const;
    void                DetachTurningPoints          ( VectorTurningPoints &TurningPoints );
    double              CalcDamage                   ( double dAmplitude, double dAvrg, double dCounts ) const;
    double              GetBKZ                       () const;
    double              GetBKZ                       ( const VectorStufen &Stufen ) const;
    double              GetBKZ                       ( const VectorStufenAmpl &Stufen ) const;
    double              GetShapeValue                () const;
    double              GetShapeValue                ( const VectorStufenAmpl &Stufen, double &dAele ) const;
    double              GetShapeValue                ( const VectorStufenAmpl &Stufen ) const;
    double              GetSumH                      () const;
    double              GetSumH                      ( const VectorStufenAmpl &Stufen ) const;


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
    }
};


void CRainflow::ClearTurningPoints()
{
    m_TurningPoints.clear();
    m_ctx.tp_cnt = 0;
}



void CRainflow::ZeroInit()
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

    m_ctx.residual_method       = m_eResiduum;

    Rainflow::clear_counts();
}



int CRainflow::EntriesCount() const
{
    unsigned count = 0;

    if( !Rainflow::rfm_non_zeros( &count ) )
    {
        ASSERT( false );
        return -1;
    }

    return (int)count;
}


// Nur fuer Matrix nach Residuenbehandlung!
// MakeSymmetric() beeinflusst nicht die Ergebnisse aus den Zaehlverfahren:
// KGUZ und BPZ (GetStufen) zaehlen stehende und haengende Hystereseaeste!

void CRainflow::MakeSymmetric() 
{
    if( !Rainfow::rfm_make_symmetric() )
    {
        ASSERT( false );
    }
}



void CRainflow::GetStufenNewR( VectorStufen &StufenOrig, VectorStufenAmpl &Stufen, 
                             double dNewR, double dM ) const
{
    t_stufe_ampl stufe;

    Stufen.clear();

    if( Rainflow::at_init( /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, dM, 
                           /*Sm_rig*/ 0.0, dNewR, /*R_pinned*/ true, /*symmetric*/ false ) )
    {
        for( int i = 0; i < (int)StufenOrig.size(); i++ ) 
        {
            double dTransformedAmpl;

            if( Rainflow::at_transform( StufenOrig[i].dRange / 2.0, StufenOrig[i].dAvrg, &dTransformedAmpl ) )
            {
                stufe.dAmpl   = dTransformedAmpl;
                stufe.dCounts = StufenOrig[i].dCounts;
            } else {
                stufe.dAmpl   = 0;
                stufe.dCounts = -1;  // TBD: Gibts eine bessere Loesung?
            }

            Stufen.push_back( stufe );
        }
    }

    if( !Rainflow::at_init( /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, m_AmplTrans_dM, 
                            /*Sm_rig*/ m_AmplTrans_dAvrg, /*R_rig*/ m_AmplTrans_dR, /*R_pinned*/ m_AmplTrans_mode = 2, /*symmetric*/ false ) )
    {
        ASSERT( false );
    }
}



void CRainflow::GetStufenNewMean( VectorStufen &StufenOrig, VectorStufenAmpl &Stufen, 
                                double dNewMean, double dM ) const
{
    t_stufe_ampl stufe;

    Stufen.clear();

    if( Rainflow::at_init( /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, dM, 
                           /*Sm_rig*/ dNewMean, /*R_rig*/ -1, /*R_pinned*/ false, /*symmetric*/ false ) )
    {
        for( int i = 0; i < (int)StufenOrig.size(); i++ ) 
        {
            double dTransformedAmpl;

            if( Rainflow::at_transform( StufenOrig[i].dRange / 2.0, StufenOrig[i].dAvrg, &dTransformedAmpl ) )
            {
                stufe.dAmpl   = dTransformedAmpl;
                stufe.dCounts = StufenOrig[i].dCounts;
            } else {
                stufe.dAmpl   = 0;
                stufe.dCounts = -1;  // TBD: Gibts eine bessere Loesung?
            }

            Stufen.push_back( stufe );
        }
    }

    if( !Rainflow::at_init( /*Sa*/ NULL, /*Sm*/ NULL, /*count*/ 0, m_AmplTrans_dM, 
                            /*Sm_rig*/ m_AmplTrans_dAvrg, /*R_rig*/ m_AmplTrans_dR, /*R_pinned*/ m_AmplTrans_mode = 2, /*symmetric*/ false ) )
    {
        ASSERT( false );
    }
}



void CRainflow::GetStufen( VectorStufen &Stufen ) const
{
    rfc_rfm_item_v buffer;

    Stufen.clear();

    if( Rainflow::rfm_get( buffer ) )
    {
        for( size_t i = 0; i < buffer.size(); i++ )
        {
            /* Es wird mit den Klassenmitten gerechnet */
            /* ( i + 0.5 ) - ( j + 0.5 ) --> j - i
             * ( i + 0.5 ) + ( j + 0.5 ) --> i + j + 1
             */
            unsigned from   = buffer[i].from;
            unsigned to     = buffer[i].to;
            double   counts = (double)buffer[i].counts / Rainflow::full_inc();
            double   range  = m_class_width * Abs( (int)from - to );
            double   avrg   = m_class_width * ( from + to + 1 ) / 2.0 + m_range_min;

            t_stufe stufe;

            stufe.dRange  = range;
            stufe.dAvrg   = avrg;
            stufe.dCounts = counts;

            Stufen.push_back( stufe );
        }
    }
    else
    {
        ASSERT( false );
    }
}



void CRainflow::GetStufen( VectorStufenAmpl &Stufen ) const
{
    rfc_counts_v rp;
    rfc_value_v  Sa;

    Stufen.clear();

    if( Rainflow::rp_get( rp, Sa ) )
    {
        for( size_t i = 0; i < rp.size(); i++ )
        {
            t_stufe_ampl stufe;

            stufe.dAmpl   = Sa[i];
            stufe.dCounts = (double)rp[i] / Rainflow::full_inc();

            Stufen.push_back( stufe );
        }
    }
    else
    {
        ASSERT( false );
    }
}



void CRainflow::GetStufen( VectorStufenFromTo &Stufen ) const
{
    rfc_rfm_item_v buffer;

    Stufen.clear();

    if( Rainflow::rfm_get( buffer ) )
    {
        for( size_t i = 0; i < buffer.size(); i++ )
        {
            double     from = m_class_width * ( 0.5 + buffer[i].from ) + m_range_min;
            double       to = m_class_width * ( 0.5 + buffer[i].to )   + m_range_min;
            double   counts = (double)buffer[i].counts / Rainflow::full_inc();
            
            t_stufe_from_to stufe;

            stufe.dFrom   = from;
            stufe.dTo     = to;
            stufe.dCounts = counts;

            Stufen.push_back( stufe );
        }
    }
    else
    {
        ASSERT( false );
    }
}



void CRainflow::GetStufen( VectorStufenKGUZ &Stufen ) const
{
    rfc_value_v  lc;
    rfc_counts_v level;

    Stufen.clear();

    if( Rainflow::lc_get( lc, level ) )
    {
        for( size_t i = 0; i < buffer.size(); i++ )
        {
            t_stufe_KGUZ stufe;

            stufe.dLevel = level[i];
            stufe.dCounts = (double)lc[i] / Rainflow::full_inc();

            Stufen.push_back( stufe );
        }
    }
    else
    {
        ASSERT( false );
    }
}



void CRainflow::SetStufen( const VectorStufen &Stufen )
{
    size_t i;

    if( !m_isParametrized || !m_rfc_ctx.rp || !m_rfc_ctx.lc ) 
    {
        ASSERT( false );
        return;
    }

    Rainflow::clear_counts();

    // Rangepair
    for( i = 0; i < Stufen.size(); i++ ) 
    {
        t_stufe stufe = Stufen[i];
        int from = ClassNo( (rfc_value_t)( stufe.dAvrg - stufe.dRange / 2.0 ) );
        int   to = ClassNo( (rfc_value_t)( stufe.dAvrg + stufe.dRange / 2.0 ) );

        ASSERT( from >= 0 && from < m_class_count 
                && to >=0 && to < m_class_count );

        AddCycles( from, to, (int)( stufe.dCounts + 0.5 ) );
    }

    if( !Rainflow::rp_from_rfm( m_rfc_ctx.rp, NULL, NULL ) ||
        !Rainflow::lc_from_rfm( m_rfc_ctx.lc, NULL, NULL, RF::RFC_FLAGS_COUNT_LC ) ||
        !Rainflow::damage_from_rfm( NULL, &m_rfc_ctx.damage ) )
    {
        ASSERT( false );
    }
}



void CRainflow::SetStufen( const VectorStufenFromTo &Stufen )
{
    size_t i;

    if( !m_isParametrized || !m_rfc_ctx.rp || !m_rfc_ctx.lc ) 
    {
        ASSERT( false );
        return;
    }

    Rainflow::clear_counts();

    for( i = 0; i < Stufen.size(); i++ ) 
    {
        t_stufe_from_to stufe = Stufen[i];
        int from = stufe.dFrom;
        int   to = stufe.dTo;

        ASSERT( from >= 0 && from < m_class_count 
                && to >=0 && to < m_class_count );

        AddCycles( from, to, (int)( stufe.dCounts + 0.5 ) );
    }

    if( !Rainflow::rp_from_rfm( m_rfc_ctx.rp, NULL, NULL ) ||
        !Rainflow::lc_from_rfm( m_rfc_ctx.lc, NULL, NULL, RF::RFC_FLAGS_COUNT_LC ) ||
        !Rainflow::damage_from_rfm( NULL, &m_rfc_ctx.damage ) )
    {
        ASSERT( false );
    }
}



void CRainflow::GetResiduum( VectorResiduum &Residuum ) const
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
}


const typename CRainflow::VectorTurningPoints& CRainflow::GetTurningPoints() const
{
    return Rainflow::rfc_tp_storage();
}



void CRainflow::GetTurningPoints( VectorTurningPoints &TurningPoints, size_t left, size_t right ) const
{
    if( !left && !right )
    {
        VectorTurningPoints aCopy( m_tp );
        TurningPoints.swap( aCopy );
    }
    else
    {
        typename VectorTurningPoints::const_iterator it;
        size_t count = 0;

        TurningPoints.clear();

#if !HAVE_NEOLIB
        if( !m_tp.empty() )
        {
            if( !left )  left  = m_tp.front().pos;
            if( !right ) right = m_tp.back().pos;
        }

        // Count number of values first...
        for( it = m_tp.begin(); it != m_tp.end(); it++ )
        {
            if( it->pos >= left && it->pos <= right )
            {
                count++;
            }
        }

        // ...then allocate space 
        TurningPoints.reserve( count );
#endif


        for( it = m_tp.begin(); it != m_tp.end(); it++ )
        {
            if( it->pos >= left && it->pos <= right )
            {
                TurningPoints.push_back( *it );
            }
        }
    }
}



void CRainflow::Parametrize( double range_min,         double range_max,
                             double range_fixed_min,   double range_fixed_max, double range_fixed, 
                             double class_count,       double class_width, 
                             double dilation,          double hysteresis ) 
{
    if( class_count > 0.0 ) 
    {
        class_count = floor( class_count + 0.5 );
    } else {
        class_count = DBL_NAN;
    }

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
    }

    /*
     *
     *   range
     *
     */
    if( !is_range_fixed ) 
    {
        range = static_cast<rfc_value_t>( class_width ) * class_count;
        
        if( !DBL_ISFINITE( range ) ) 
        {
            if( is_range_max_fixed && is_range_min_fixed ) {
                range = Abs( static_cast<rfc_value_t>( range_fixed_max ) - static_cast<rfc_value_t>( range_fixed_min ) );
            } else if( is_range_max_fixed && !is_range_min_fixed ) {
                range = Abs( static_cast<rfc_value_t>( range_fixed_max ) - static_cast<rfc_value_t>( range_min ) );
            } else if( !is_range_max_fixed && is_range_min_fixed ) {
                range = Abs( static_cast<rfc_value_t>( range_max ) - static_cast<rfc_value_t>( range_fixed_min ) );
            } else {
                range = Abs( static_cast<rfc_value_t>( range_max ) - static_cast<rfc_value_t>( range_min ) );
            }
        }
    } else {
        range = range_fixed;
    }

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
            range_min = static_cast<rfc_value_t>( range_fixed_max ) - static_cast<rfc_value_t>( range );
        } else if( DBL_ISFINITE( range_max ) ) {
            range_min = static_cast<rfc_value_t>( range_max ) - static_cast<rfc_value_t>( range );
        }
    } else {
        range_min = range_fixed_min;
    }
    
    if( !is_range_max_fixed ) 
    {
        if( is_range_min_fixed ) 
        {
            range_max = static_cast<rfc_value_t>( range_fixed_min ) + static_cast<rfc_value_t>( range );
        } else if( DBL_ISFINITE( range_min ) ) {
            range_max = static_cast<rfc_value_t>( range_min ) + static_cast<rfc_value_t>( range );
        }
    } else {
        range_max = range_fixed_max;
    }

    ASSERT( DBL_ISFINITE ( range_min + range_max ) );
    ASSERT( range_max >= range_min );
    range = range_max - range_min;

    if( is_range_max_fixed && is_range_min_fixed ) 
    {
        is_range_fixed = true;
    }
    
    /*
     *
     *   class_count
     *
     */
    if( is_class_width_fixed && !is_class_count_fixed ) 
    {
        ASSERT( class_width > 0.0 );
        class_count = (int)Ceil( static_cast<rfc_value_t>( range ) / static_cast<rfc_value_t>( class_width ) );
    } /* end if */
    
    if( !is_class_count_fixed && !is_class_width_fixed ) 
    {
        class_count = DEFAULT_CLASS_COUNT;
        is_class_count_fixed = 1;
    }

    if( !is_dilation_fixed ) 
    {
        dilation = (double)DEFAULT_DILATION / 100.0;
        is_dilation_fixed = 1;
    }
    
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
    }
    
    
    /* Einfache Klassenbreite, 2x halbe Klassenbreite Zuschlag, wenn dilation_trunc==1 */
    if( !is_class_width_fixed ) class_width = static_cast<rfc_value_t>( range ) / IROUND( class_count - dilation_trunc );
    /* Halbe Klassenbreite auch an den Klassiergrenzen aufschlagen */
    if( !is_range_min_fixed )   range_min  -= static_cast<rfc_value_t>( class_width ) * dilation_trunc / 2.0f;
    if( !is_range_max_fixed )   range_max  += static_cast<rfc_value_t>( class_width ) * dilation_trunc / 2.0f;
    /* Klassierbereich neu bestimmen */
    if( !is_range_fixed)        range       = static_cast<rfc_value_t>( range_max ) - static_cast<rfc_value_t>( range_min );

    // Dilation unabhaengig von fixierten Werten immer anwenden (muss auf 0 gesetzt werden, wenn nicht gewuenscht!)
    range_min -= static_cast<rfc_value_t>( range )     * static_cast<rfc_value_t>( dilation_frac ) / 2.0f;
    range_max += static_cast<rfc_value_t>( range )     * static_cast<rfc_value_t>( dilation_frac ) / 2.0f;
    range      = static_cast<rfc_value_t>( range_max ) - static_cast<rfc_value_t>( range_min );

    
    if( DEFAULT_ROUNDOFF > 0 ) 
    {
        range_min = static_cast<rfc_value_t>( floor( static_cast<rfc_value_t>( range_min ) * DEFAULT_ROUNDOFF ) / DEFAULT_ROUNDOFF );
        range_max = static_cast<rfc_value_t>( floor( static_cast<rfc_value_t>( range_max ) * DEFAULT_ROUNDOFF ) / DEFAULT_ROUNDOFF );
        range = static_cast<rfc_value_t>( range_max ) - static_cast<rfc_value_t>( range_min );
    }

    if( !is_class_width_fixed ) class_width = static_cast<rfc_value_t>( range ) / IROUND( class_count );
    if( !is_class_count_fixed ) class_count = (int)Ceil( static_cast<rfc_value_t>( range ) / static_cast<rfc_value_t>( class_width ) );

    if( !DBL_ISFINITE( hysteresis ) || hysteresis == DEFAULT_HYSTERESIS ) 
    {
        hysteresis = class_width;
    }

    ASSERT( DBL_ISFINITE( hysteresis ) && hysteresis >= 0.0 );

    m_class_count = IROUND( class_count );
    m_class_width = static_cast<rfc_value_t>( class_width );
    m_hysteresis  = static_cast<rfc_value_t>( hysteresis );
    m_range_min   = static_cast<rfc_value_t>( range_min );
    m_dilation    = static_cast<rfc_value_t>( dilation );

    if( m_class_width <= 0.0 )
    {
        m_class_width = 10.0 * DBL_EPSILON;
    }

    m_isParametrized = ( m_class_count > 0 ) && ( m_class_width > 0.0 );

    if( m_isParametrized )
    {
        Rainflow::init( m_class_count, m_class_width, m_range_min, m_hysteresis, RFC_FLAGS_DEFAULT );

        SetSNCurve( m_SNCurve.SD, m_SNCurve.ND, m_SNCurve.k, m_SNCurve.k2 );

        switch( m_AmplTrans_mode )
        {
            case 0:
                SetAmplTrans( DBL_NAN, DBL_NAN, false );
                break;
            case 1:
                SetAmplTrans( m_AmplTrans_dM, m_AmplTrans_dR, true );
                break;
            case 2:
                SetAmplTrans( m_AmplTrans_dM, m_AmplTrans_dAvrg, false );
                break;
            default:
                ASSERT( false );
        }
    }
}



void CRainflow::DoRainflow( const double *stream_in, size_t in_len, bool do_finalize ) 
{
    const   double*             pin;                        /* Pointer auf Zeitreihe (Eingangsdaten) */
            double              dummy = 0.0;
            size_t              pos;                        /* Aktuelle Position in der Zeitreihe (pin) */
            rfc_value_t                   d0, dl;                     /* Rueckwaertsgerichtete Gradienten ueber die letzten 3 Punkte (res[end-1],res[end],*pin) */

    if( !m_isParametrized || ( !stream_in && in_len ) ) 
    {
        ASSERT( false );
        return;
    }

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
    Rainflow::feed( pin, in_len );

    if( do_finalize )
    {
        switch( m_eResiduum )
        {
            case RESIDUUM_NONE:
            case RESIDUUM_IGNORE:
            case RESIDUUM_HALFCYCLES:
            case RESIDUUM_FULLCYCLES:
            case RESIDUUM_REPEATED:
            case RESIDUUM_CLORMANN_SEEGER:
            case RESIDUUM_RP_DIN:
                Rainflow::finalize( m_eResiduum ); break;
            default:
                ASSERT( false );
                return;
        }
    }
}



void CRainflow::DoRainflow( const VectorData& stream ) 
{
    if( !stream.empty () ) 
    {
        DoRainflow( &stream[0], stream.size() );
    }
    else
    {
        double dDummy = 0.0;

        DoRainflow( &dDummy, /*in_len*/ 0, true );
    }
}



template<class InputIt>
void CRainflow::DoRainflow( InputIt first, InputIt last )
{
    size_t distance = 0;
    std::vector<rfc_value_t> stride( 1024 );
    
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


 
double CRainflow::CalcDamage( double dAmplitude, double dAvrg, double dCounts ) const
{
    double dDamage;
    double dLW;

    Rainflow::at_transform( dAmplitude, dAvrg, &dAmplitude );

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
}



double CRainflow::GetBKZ( const VectorStufen& Stufen ) const 
{
    double dDamageSum = 0.0;

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    }

    for( int i = 0; i < (int)Stufen.size(); i++ ) 
    {
        dDamageSum += CalcDamage( Stufen[i].dRange / 2.0, Stufen[i].dAvrg, Stufen[i].dCounts );
    }

    return dDamageSum;
}



double CRainflow::GetBKZ( const VectorStufenAmpl& Stufen ) const 
{
    double dDamageSum = 0.0;

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    }

    for( int i = 0; i < (int)Stufen.size(); i++ ) 
    {
        dDamageSum += CalcDamage( Stufen[i].dAmpl, 0.0 /* dAvrg */, Stufen[i].dCounts );
    }

    return dDamageSum;
}



double CRainflow::GetBKZ() const 
{
    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    }

    VectorStufen Stufen;
    GetStufen( Stufen );

    return GetBKZ( Stufen );
}



double CRainflow::GetSumH( const VectorStufenAmpl& Stufen ) const 
{
    double dAmount = 0.0;

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    }

    for( int i = 0; i < (int)Stufen.size(); i++ ) 
    {
        dAmount += Stufen[i].dCounts;
    }

    return dAmount;
}



double CRainflow::GetSumH() const 
{
    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    } 

    VectorStufenAmpl Stufen;
    GetStufen( Stufen );

    return GetSumH( Stufen );
}



double CRainflow::GetShapeValue( const VectorStufenAmpl& Stufen, double& dAele ) const 
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
    }

    // Maximale Amplitude im Kollektiv suchen
    for( int i = 0; i < (int)Stufen.size(); i++ ) 
    {
        if( Stufen[i].dCounts > 0 && Stufen[i].dAmpl > dMaxAmpl )
        {
            dMaxAmpl = Stufen[i].dAmpl;
        }
    }

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
    }

    dAele = 1.0 / dV;  // Abstand zwischen Bauteil- und Lebensdauerlinie, Miner-elementar

    // Voelligkeitsmass zurueckgeben
    return pow( dV, 1.0 / dSNk1 );  // TODO
}



double CRainflow::GetShapeValue( const VectorStufenAmpl& Stufen ) const 
{
    double dAele;

    return GetShapeValue( Stufen, dAele );
}



double CRainflow::GetShapeValue() const 
{
    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return -1.0;
    }

    VectorStufenAmpl Stufen;
    GetStufen( Stufen );
    
    return GetShapeValue( Stufen );
}
