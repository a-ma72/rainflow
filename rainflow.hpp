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

#include <vector>
#include <algorithm>
#include "rainflow.hpp"

#undef ASSERT
#include <assert.h>
#define ASSERT(x) assert(x)

using namespace rainflowC as RF;

template<class T> class CRainflowT;
typedef CRainflowT<RF::RFC_value_type> CRainflow;

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
        DEFAULT_HYSTERESIS  = -1,        /* Defaultwert wird zu jew. Klassenbreite! */
        DEFAULT_CLASS_COUNT = 100,       /* Default fuer EGDB Streckenvergleich */
        DEFAULT_DILATION    = 110,       /* Default fuer EGDB Streckenvergleich (in Prozent!)*/
        DEFAULT_WL_SD       = 1000,      /* Default fuer EGDB Streckenvergleich --> 1E3 */
        DEFAULT_WL_ND       = 10000000,  /* Default fuer EGDB Streckenvergleich --> 1E7 */
        DEFAULT_WL_k        = -5,        /* Default fuer EGDB Streckenvergleich */
        DEFAULT_ROUNDOFF    = 1000000    /* Fuer Vergleichbarkeit mit LMS TecWare */
    }; // end enum

    static const int HALF_CYCLE_INCREMENT = 1;
    static const int FULL_CYCLE_INCREMENT = 2;

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


public:
    explicit 
    CRainflowT()
    : m_iProgressState( -1 )
    { 
        RF::rfc_ctx_s dummy = { sizeof RF::rfc_ctx_s };

        m_SNCurve.SD = DEFAULT_WL_SD;
        m_SNCurve.ND = DEFAULT_WL_ND;
        m_SNCurve.k1 = DEFAULT_WL_k;
        m_SNCurve.k2 = DEFAULT_WL_k;

        m_rfc_ctx = dummy;

        m_doSpreadDamage = SPRDAM_HALF_23;  

        SetAmplTrans( DBL_NAN, DBL_NAN, false );
        ZeroInit(); 
    } // end of CRainflowT
    

    virtual ~CRainflowT()                                            
    { 
        RF::RFC_deinit( &m_rfc_ctx );
    } // end of ~CRainflowT                                         
    

    // Gibt Anzahl Schwingspiele x FULL_CYCLE_INCREMENT zurueck!
    int GetCountIncrements( int from, int to ) const
    {
        RF::RFC_counts_type counts;

        ASSERT( m_isParametrized );

        if( RF::RFC_rfm_peek( &m_rfc_ctx, from, to, &counts ) )
        {
            return counts;
        }

        return -1;
    } // end of GetCountIncrements
    

    void SetCounts( int from, int to, int counts )
    {
        RF::RFC_counts_type counts_;

        ASSERT( counts >= 0 );
        ASSERT( m_isParametrized );
        ASSERT( from < m_class_count && from >= 0 && to < m_class_count && to >= 0 );

        counts_ = (RF::RFC_counts_type)counts;

        if( m_isSymmetric && to < from ) 
        {
            int temp = from;
            to = from;
            from = temp;
        } // end if

        RF::RFC_rfm_poke( &m_rfc_ctx, from, to, counts_ * FULL_CYCLE_INCREMENT );
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

        RF::RFC_wl_init_original( &m_rfc_ctx, SD, ND, k );
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

        if( m_isSymmetric && to < from ) 
        {
            int temp = from;
            to = from;
            from = temp;
        } // end if

        RF::RFC_rfm_poke( &m_rfc_ctx, from, to, cycles * FULL_CYCLE_INCREMENT, /*add_only*/ true );
    } // end of AddCycles
    

    inline 
    void AddHalfCycle( int from, int to )
    {
        ASSERT( m_isParametrized );
        ASSERT( from < m_class_count && from >= 0 && to < m_class_count && to >= 0 );

        if( m_isSymmetric && to < from ) 
        {
            int temp = from;
            to = from;
            from = temp;
        } // end if

        RF::RFC_rfm_poke( &m_rfc_ctx, from, to, cycles * HALF_CYCLE_INCREMENT, /*add_only*/ true );
    } // end of AddHalfCycle
    
    
    double    NoValue                      () const { return DBL_NAN; }
    void      ClearResiduum                ();
    void      ClearTurningPoints           ();
    void      FreeMatrix                   ( bool bFreeAll = true );
    void      ZeroInit                     ();
    void      PushMatrix                   ( t_matrix_backup *pMatrixBackup = NULL ); 
    void      RestoreMatrix                ( t_matrix_backup *pMatrixBackup = NULL );
    void      ReDimMatrix                  ();
    void      EntriesCount                 () const;
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
void CRainflowT<T>::FreeMatrix( bool bFreeAll ) 
{
    if ( bFreeAll && m_rf_matrix_no_res ) 
    {
        free( m_rf_matrix_no_res );
        m_rf_matrix_no_res = NULL;
    } /* end if */

    if ( m_rf_matrix ) 
    {
        free( m_rf_matrix );
        m_rf_matrix = NULL;
        m_eResiduum = RESIDUUM_NONE;
    } /* end if */
} // end of FreeMatrix


template<class T>
void CRainflowT<T>::ClearResiduum( int iStack /*=-1*/ )
{
    if( iStack < 0 )
    {
        m_Residuum[0].clear();
        m_Residuum[1].clear();
    } else {
        m_Residuum[iStack].clear();
    }
} // end of ClearResiduum


template<class T>
void CRainflowT<T>::ClearTurningPoints( int iStack /*=-1*/ )
{
    if( iStack < 0 )
    {
        m_TurningPoints[0].clear();
        m_TurningPoints[1].clear();
    } else {
        m_TurningPoints[iStack].clear();
    }
} // end of ClearTurningPoints


template<class T>
void CRainflowT<T>::ZeroInit()
{
    RestoreMatrix();
    FreeMatrix();
    ClearResiduum();
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
    m_min                       =  EmptyRes();
    m_max                       =  EmptyRes();
} // end of ZeroInit


template<class T>
void CRainflowT<T>::PushMatrix( t_matrix_backup *pMatrixBackup ) 
{ 
    ASSERT( m_isParametrized );

    if( pMatrixBackup == NULL ) {
        pMatrixBackup = &m_matrix_backup;
    } // end if

    t_matrix_backup backup = { m_rf_matrix_no_res, m_rf_matrix, m_eResiduum };

    *pMatrixBackup = backup;

    m_rf_matrix = NULL; 
    m_rf_matrix_no_res = NULL; 
    m_eResiduum = RESIDUUM_NONE;
} // end of PushMatrix


template<class T>
void CRainflowT<T>::RestoreMatrix( t_matrix_backup *pMatrixBackup )
{
    if( pMatrixBackup == NULL ) 
    {
        pMatrixBackup = &m_matrix_backup;
    } // end if

    t_matrix_backup null = { 0 };

    FreeMatrix();
    m_rf_matrix_no_res = pMatrixBackup->m_rf_matrix_no_res;
    m_rf_matrix = pMatrixBackup->m_rf_matrix;
    m_eResiduum = pMatrixBackup->m_eResiduum;

    *pMatrixBackup = null;
} // end of RestoreMatrix


template<class T>
void CRainflowT<T>::ReDimMatrix()
{
    ASSERT( m_isParametrized );

    FreeMatrix();
    
    m_rf_matrix_no_res = (int*)calloc( m_class_count * m_class_count, sizeof( int ) );
} // end of ReDimMatrix


template<class T>
void CRainflowT<T>::EntriesCount() const
{
    int i, j;
    int entries_count = 0;

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
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
    int i, j;
    int *prf = GetMatrix();

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
    } // end if

    ASSERT( prf );
    
    /* Matrix entlang der Diagonalen "falten" */
    for ( i = 0; i < m_class_count; i++ )   // from
    {
        for ( j = i + 1; j < m_class_count; j++ )   // to
        {
            // j immer > i !
            prf[ i * m_class_count + j ] += prf[ j * m_class_count + i ];
            prf[ j * m_class_count + i ] = 0;  // Negative Aeste zu Null setzen
        } /* end for */
    } /* end for */
} /* end of MakeSymmetric() */


template<class T>
void CRainflowT<T>::GetStufenNewR( VectorStufen &StufenOrig, VectorStufenAmpl &Stufen, 
                                         double dNewR, double dM ) const
{
    AmplTransform AmplTransform;
    t_stufe_ampl stufe;

    Stufen.clear();

    for( int i = 0; i < (int)StufenOrig.size(); i++ ) 
    {
        double dTransformedAmpl;

        dTransformedAmpl = AmplTransform.Calc( StufenOrig[i].dRange / 2.0, 
                                               StufenOrig[i].dAvrg, 
                                               dM, dNewR, true );
        if( DBL_ISFINITE( dTransformedAmpl ) )
        {
            stufe.dAmpl   = dTransformedAmpl;
            stufe.dCounts = StufenOrig[i].dCounts;
        } else {
            stufe.dAmpl   = 0;
            stufe.dCounts = -1;  // TBD: Gibts eine bessere Loesung?
        } // end if

        Stufen.push_back( stufe );
    } // end for
} // end of GetStufenNewR


template<class T>
void CRainflowT<T>::GetStufenNewMean( VectorStufen &StufenOrig, VectorStufenAmpl &Stufen, 
                                            double dNewMean, double dM ) const
{
    AmplTransform AmplTransform;
    t_stufe_ampl stufe;

    Stufen.clear();

    for( int i = 0; i < (int)StufenOrig.size(); i++ ) 
    {
        double dTransformedAmpl;

        dTransformedAmpl = AmplTransform.Calc( StufenOrig[i].dRange / 2.0, 
                                               StufenOrig[i].dAvrg, 
                                               dM, dNewMean, false );
        if( DBL_ISFINITE( dTransformedAmpl ) )
        {
            stufe.dAmpl   = dTransformedAmpl;
            stufe.dCounts = StufenOrig[i].dCounts;
        } else {
            stufe.dAmpl   = 0;
            stufe.dCounts = -1;  // TBD: Gibts eine bessere Loesung?
        } // end if

        Stufen.push_back( stufe );
    } // end for
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

        if ( counts > 0.0 ) 
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

    if( m_rf_matrix_no_res ) 
    {
        for( int i = 1; i < m_class_count; i++ ) 
        {
            double dCounts = 0.0;
            
            // Erster Index [0] zaehlt das Durchschreiten der ersten oberen Klassengrenze.
            KGUZ_Temp[ i - 1 ].dLevel = i * m_class_width + m_range_min;
            
            for( int j = i; j < m_class_count; j++ )  // To
            {
                for( int k = 0; k < i; k++ )          // From
                {
                    // Immer symmetrisch berechnen (Schleife deckt nur die Haelfte (Dreieck) ab!)
                    // Wichtiger Hinweis: 
                    // Die hier verwendete Matrix ist keine reine(!) Uebergangsmatrix wie z.B. in der 
                    // FVA Richtlinie "Zaehlverfahren 2010" beschrieben! Da geschlossene Hysteresen gezaehlt werden
                    // gibt es zu einer negativen Flanke auch immer eine positive. Daher muessen die Haeufigkeiten
                    // aus beiden Matrixhaelften beruecksichtigt werden:
                    dCounts += (double) m_rf_matrix_no_res[ k * m_class_count + j ] / FULL_CYCLE_INCREMENT;
                    dCounts += (double) m_rf_matrix_no_res[ j * m_class_count + k ] / FULL_CYCLE_INCREMENT;
                } // end for
            } // end for

            KGUZ_Temp[ i - 1 ].dCounts = dCounts;
        } // end for
    } // end if

    /* Residuum einrechnen */
    if( m_Residuum[0].size() >= 2 ) 
    {
        t_res from, to;

        from = m_Residuum[0][0];
        for ( size_t n = 1; n < m_Residuum[0].size(); n++ ) 
        {
            to = m_Residuum[0][n];
            if( from.iCno < to.iCno ) 
            {
                /* Nur die aufsteigenden Aeste beruecksichtigen */
                for( int i = from.iCno; i < to.iCno; i++ ) 
                {
                    KGUZ_Temp[ i ].dCounts += 1.0;
                } /* end for */
            } /* end if */
            from = to;
        } /* end for */
    } /* end if */

    for ( size_t n = 0; n < KGUZ_Temp.size(); n++ ) 
    {
        if( KGUZ_Temp[ n ].dCounts > 0.0 ) 
        {
            Stufen.push_back( KGUZ_Temp[ n ] );
        } // end if
    } // end for
} // end of GetStufen


template<class T>
void CRainflowT<T>::SetStufen( const VectorStufen &Stufen )
{
    size_t i;

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
    } // end if

    ReDimMatrix();

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
} // end of SetStufen


template<class T>
void CRainflowT<T>::SetStufen( const VectorStufenFromTo &Stufen )
{
    size_t i;

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
    } // end if

    ReDimMatrix();

    for( i = 0; i < Stufen.size(); i++ ) 
    {
        t_stufe_from_to stufe = Stufen[i];
        int from = stufe.dFrom;
        int   to = stufe.dTo;

        ASSERT( from >= 0 && from < m_class_count 
                && to >=0 && to < m_class_count );

        AddCycles( from, to, (int)( stufe.dCounts + 0.5 ) );
    } // end for
} // end of SetStufen


template<class T>
void CRainflowT<T>::GetResiduum( VectorResiduum &Residuum ) const
{
    Residuum = m_Residuum[0];
} // end of GetResiduum


template<class T>
void CRainflowT<T>::SetResiduum( const VectorResiduum &Residuum )
{
    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
    } // end if

    m_Residuum[0] = Residuum;
} // end of SetResiduum


template<class T>
void CRainflowT<T>::GetResiduumStufen( VectorStufen &Stufen, e_residuum how ) const
{
    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
    } // end if

    int *backup[] = { m_rf_matrix_no_res, m_rf_matrix };
    e_residuum old_residuum = m_eResiduum;

    ReDimMatrix();
    CountResiduum( how );
    GetStufen( Stufen );
    FreeMatrix();

    m_eResiduum = old_residuum;
    m_rf_matrix_no_res = backup[0];
    m_rf_matrix = backup[1];
} // end of GetResiduumStufen


template<class T>
void CRainflowT<T>::GetResiduumStufen( VectorStufenFromTo &Stufen, e_residuum how ) const
{
    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
    } // end if

    int *backup[] = { m_rf_matrix_no_res, m_rf_matrix };
    e_residuum old_residuum = m_eResiduum;

    ReDimMatrix();
    CountResiduum( how );
    GetStufenFromTo( Stufen );
    FreeMatrix();

    m_eResiduum = old_residuum;
    m_rf_matrix_no_res = backup[0];
    m_rf_matrix = backup[1];
} // end of GetResiduumStufen


template<class T>
const typename CRainflowT<T>::VectorTurningPoints& CRainflowT<T>::GetTurningPoints() const
{
    return m_TurningPoints[0];
}


template<class T>
void CRainflowT<T>::DetachTurningPoints( VectorTurningPoints &TurningPoints )
{
    TurningPoints.swap( m_TurningPoints[0] );
    m_TurningPoints[0].clear();
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


template<class T>
void CRainflowT<T>::Parametrize( double range_min,         double range_max,
                                       double range_fixed_min,   double range_fixed_max, double range_fixed, 
                                       double class_count,       double class_width, 
                                       double dilation,          double hysteresis ) 
{
    m_min = m_max = EmptyRes();

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
            VectorResiduum&     Residuum( m_Residuum[0] );  /* Residuum */

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
        m_min = m_max = EmptyRes();

        ReDimMatrix();             /* Neue Matrix anlegen */
        ClearResiduum(0);          /* Altes Residuum loeschen */
        ClearTurningPoints(0);     /* Umkehrpunkte loeschen */
        
        m_bCountPending = !do_finalize;
    }
    else
    {
        if( Residuum.size() > 1 )
        {
            // Gradient am Ende des Residuums
            d0 = Residuum.back().Value - ( Residuum.rbegin() + 1 )->Value;
        } // end if
    } // end if
    
    ASSERT( m_rf_matrix_no_res && pin );

    /* Schleife ueber alle Punkte */
    for( pos = 0; pos <= in_len; pos++, pin++ ) 
    {
        /* Letzter Messpunkt im Arbeitsstack ist kein endgueltiger Umkehrpunkt.
        /* Da der Aufruf mehrfach fuer kontinuierliche Messsignale
        /* aufgerufen werden kann, muss kenntlich gemacht werden, wenn
        /* keine weiteren Daten reingereicht werden. In diesem Fall muss
        /* der letzte Punkt des Arbeitsstackes als endgueltiger Umkehrpunkt
        /* betrachtet werden (do_finalize).
        */
        bool    read_pin = true;
        size_t  nIdx     = m_CurrentPos + pos + 1;  /* Index, base 1 */
        size_t  tp_last  = 0;
        
        if( pos == in_len ) 
        {
            read_pin = false;
            if( !do_finalize ) 
            {
                /* Wenn noch weitere Daten folgen werden, dann hier abbrechen */
                break;
            } /* end if */
        } /* end if */

        if( !read_pin ) 
        {
            size_t st_count = Residuum.size();

            /* lokales Ende des Streams erreicht, tp_last == 0 */
            if( do_finalize && st_count )
            {
                /* Hier ist das (globale) Ende des Streams erreicht. Der Letzte Punkt wird noch aktualisiert und
                   als Umkehrpunkt eingetragen.
                   Aber nur aktualisieren, wenn sich der Punkt von dem voherigen unterscheidet */
                T value = Residuum.back().Value;
                Residuum.back().iCno = ClassNo( value );

                t_turning_point tp;
                tp.Value   = Residuum.back().Value;
                tp.Avrg    = tp.Value;
                tp.nIdx    = Residuum.back().nIdx;  /* Index, base 1 */
                tp.nIdxAdj = 0;                     /* Index, base 1 */
                tp.BKZ     = 0.0; 
                AddTurningPoint( tp );  // Veraendert Residuum.back()

                tp_last = Residuum.size(); /* Letzter verwertbarer Umkehrpunkt im Stack, base 1 */
            }
        } 
        else 
        {  /* if( !read_pin ) */
            size_t st_count = Residuum.size();
            
            if( !DBL_ISFINITE(*pin) ) 
            {
                continue;
            }

            /* Naechsten Messpunkt holen.
             * Mind. 2 Punkte werden fuer die Suche nach turning points benoetigt.
             * Wenn der Arbeitsstack 2 Punkte aufgenommen hat, werden es nie wieder
             * weniger als 2 sein, da eine geschlossene Hysterese nur durch mind. 4
             * Datenpunkte erkannt werden kann und dann auch nur 2 Punkte entfernt werden.
             */
            if( st_count < 2 ) 
            {
                t_res res;

                res.Value   = static_cast<T>( *pin );
                res.iCno    = ClassNo( res.Value );
                res.nIdx    = nIdx;
                res.nIdx_TP = 0;                       /* Index Turning Points, base 1 */

                /* 
                 * Noch kein Umkehrpunkt im Stack
                 */
                if( !m_min.nIdx )
                {
                    m_min = m_max = res;
                }
                else
                {
                    if( res.Value < m_min.Value )
                    {
                        m_min = res;
                    }
                    else if( res.Value > m_max.Value )
                    {
                        m_max = res;
                    }

                    if( m_max.Value - m_min.Value > m_hysteresis )
                    {
                        t_turning_point tp;

                        tp.nIdxAdj = 0;            /* Index, base 1 */
                        tp.BKZ     = 0.0; 

                        if( m_min.nIdx < m_max.nIdx )
                        {
                            Residuum.push_back( m_min );
                            tp.Value = m_min.Value;
                            tp.Avrg  = m_min.Value;
                            tp.nIdx  = m_min.nIdx;     /* Index, base 1 */
                            AddTurningPoint( tp );     /* Veraendert Residuum.back() */
                            Residuum.push_back( m_max );
                        }
                        else
                        {
                            Residuum.push_back( m_max );
                            tp.Value = m_max.Value;
                            tp.Avrg  = m_max.Value;
                            tp.nIdx  = m_max.nIdx;     /* Index, base 1 */
                            AddTurningPoint( tp );     /* Veraendert Residuum.back() */
                            Residuum.push_back( m_min );
                        }
                        // Gradient am Ende des Residuums
                        d0 = Residuum.back().Value - ( Residuum.rbegin() + 1 )->Value;
                    }
                }
        
                continue;
            } /* end if */

            /* Ab hier sind mind. 2 Punkte im Arbeitsstack */
            dl = (double)(static_cast<T>( *pin ) - Residuum.back().Value);  /* Delta zum letzten Umkehrpunkt */

            if ( Sign( dl ) == Sign( d0 ) ) 
            {
                /*
                Der urspruengliche Verlauf wird fortgesetzt.
                Amplitude ist automatisch groesser und braucht daher nicht ueberprueft
                zu werden!
                */

                /*
                 *
                 *        *     pin[0]
                 *       /      dl
                 *      *       stack[st_count-1]
                 *   \ /        d0
                 *    *         stack[st_count-2]
                 *
                 */

                /* Letzten lokalen Umkehrpunkt korrigieren */
                Residuum.back().nIdx   = nIdx; /* Index, base 1 */
                Residuum.back().Value  = static_cast<T>( *pin );
                Residuum.back().iCno   = -1; /* speed improvement reason */

                /* Es hat sich nur der letzte Punkt im Arbeitsstack geaendert
                /* weiter geht es mit dem naechsten Messpunkt...
                 */
                continue;

            } else {
                int cno_this;

                if( Abs( dl ) <= m_hysteresis ) 
                {
                    /* Abstand zum letzten Punkt noch zu gering --> Kein neuer Umkehrpunkt */
                    continue;
                } /* end if */

                d0 = dl;

                cno_this = ClassNo( static_cast<T>( *pin ) );

                if( Residuum.back().iCno < 0 ) 
                {
                    T value = Residuum.back().Value;
                    Residuum.back().iCno = ClassNo( value );
                } /* end if */

                /*   Neuen Umkehrpunkt gefunden!
                 *
                 *                         *      stack[st_count-1]
                 *                        / \
                 *   d0                \ /   \    dl
                 *   stack[st_count-2]  *     *   pin[0]
                 *
                 */

                /* Der letzte Umkehrpunkt im Arbeitsstack wird zum Umkehrpunkt */
                t_turning_point tp;
                tp.Value   = Residuum.back().Value;
                tp.Avrg    = tp.Value;
                tp.nIdx    = Residuum.back().nIdx;  /* Index, base 1 */
                tp.nIdxAdj = 0;                     /* Index, base 1 */
                tp.BKZ     = 0.0; 
                AddTurningPoint( tp );  // Veraendert Residuum.back()
                
                /* Neuer Wert wird neuer Umkehrpunkt im Arbeitsstack */
                t_res res;
                res.Value   = static_cast<T>( *pin );
                res.iCno    = cno_this;
                res.nIdx    = nIdx;                    /* Index, base 1 */
                res.nIdx_TP = 0;                       /* Index Turning Points, base 1 */
                Residuum.push_back( res );

                d0 = dl;

                tp_last = Residuum.size() - 1; /* Letzter verwertbarer Umkehrpunkt im Stack, base 1 (letzter Wert gehoert noch nicht dazu!)*/
            } /* end if ( Sign( dl ) == Sign( d0 ) ) */
        } /* end if ( read_pin ) */

        /* Stack hat neuen Umkehrpunkt. Daher 4P-Algorithmus
        /* zur Suche nach geschlossenen Hysteresen anstossen.
         */

        /* 4 Punkte im Arbeitsstack erforderlich -> 4P-Algorithmus */
        while( tp_last >= 4 ) 
        {
            int p1, p2, p3, p4;

            p1  = Residuum[ tp_last - 4 ].iCno;
            p2  = Residuum[ tp_last - 3 ].iCno;
            p3  = Residuum[ tp_last - 2 ].iCno;
            p4  = Residuum[ tp_last - 1 ].iCno;

            ASSERT( p1 >= 0 && p1 < m_class_count &&
                    p2 >= 0 && p2 < m_class_count &&
                    p3 >= 0 && p3 < m_class_count &&
                    p4 >= 0 && p4 < m_class_count );

            if ( Min( p2, p3 ) < Min( p1, p4 ) ||
                 Max( p2, p3 ) > Max( p1, p4 ) ) 
            {

                break;
            } /* end if */

            /* Ein Schwingspiel in die Matrix From->To eintragen */
            AddCycles( p2, p3, 1 );

            SetAdjacentTurningPoint( Residuum[ tp_last - 3 ].nIdx_TP, Residuum[ tp_last - 2 ].nIdx_TP );
            
            //int idx;
            double dAmplitude  = Abs( p2 - p3 )     * m_class_width / 2.0;
            double dAvrg       = Abs( p2 + p3 + 1 ) * m_class_width / 2.0 + m_range_min;
            
            /* Schaedigung den verantwortlichen Umkehrpunkten zuordnen */
            if( m_doSpreadDamage == SPRDAM_RAMP_DAMAGE_24 ||
                m_doSpreadDamage == SPRDAM_RAMP_AMPLITUDE_24 )
            {
                YieldDamageOverTurningPoints( Residuum[ tp_last - 3 ].nIdx_TP,   // Von P2
                                              Residuum[ tp_last - 1 ].nIdx_TP,   // Zu P4
                                              dAmplitude, dAvrg );
            } else {
                YieldDamageOverTurningPoints( Residuum[ tp_last - 3 ].nIdx_TP,   // Von P2
                                              Residuum[ tp_last - 2 ].nIdx_TP,   // Zu P3
                                              dAmplitude, dAvrg );
            }
            

            /* Umkehrpunkte p2 und p3 aus dem Residuum loeschen */
            if(0) 
            {
              for( size_t i = tp_last; i <= Residuum.size(); i++ ) 
              {
                  Residuum[ i - 3 ] = Residuum[ i - 1 ];
              } /* end for */

              Residuum.pop_back();
              Residuum.pop_back();

              tp_last -= 2;
            } else {
              Residuum.erase( Residuum.begin() + tp_last - 3, Residuum.begin() + tp_last - 1 );
              tp_last -= 2;
            } // end if

        } /* end while, 4P Alg.*/
    } /* end for, pin */

    /* Liegt noch ein Residuum vor? */
    if ( Residuum.size() ) 
    {
        typename VectorResiduum::iterator It;
        
        /* Klassenmitten nachtraeglich berechnen */
        for ( It = Residuum.begin(); It != Residuum.end(); It++ ) 
        {
            It->ClassMean = (T)( (double)It->iCno + 0.5 ) * m_class_width + m_range_min;
        } /* end for */
    } /* end if */

    m_CurrentPos += MIN( pos, in_len );

    if( do_finalize ) 
    {
        m_bCountPending = false;
        m_iProgressState = 100;
    } // end if
} // end of DoRainflow


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
void CRainflowT<T>::CountResiduum( e_residuum how ) 
{
    VectorResiduum& Residuum( m_Residuum[0] );

    if( !m_isParametrized ) 
    {
        ASSERT( false );
        return;
    } // end if

    ASSERT( m_rf_matrix_no_res );

    FreeMatrix( false ); // Nur Matrix mit Residuum loeschen
    m_eResiduum = how;

    if( how == RESIDUUM_NONE ) 
    {
        return;
    } // end if

    m_rf_matrix = (int*)calloc( m_class_count * m_class_count, sizeof( int ) );
    memcpy( m_rf_matrix, m_rf_matrix_no_res, sizeof( int ) * m_class_count * m_class_count );

    if( !Residuum.size() || how == RESIDUUM_IGNORE ) 
    {
        return;
    } // end if

    switch( m_eResiduum ) 
    {
        case RESIDUUM_HALFCYCLES:
            /* ASTM method */
            if( Residuum.size() >= 2 ) 
            {
                typename VectorResiduum::const_iterator It;
                t_res from, to;
                
                from = Residuum[0];
                for ( It = Residuum.begin() + 1; It != Residuum.end(); It++ ) 
                {
                    to = *It;
                    AddHalfCycle( from.iCno, to.iCno );
                    
                    double dAmplitude = Abs( from.iCno - to.iCno )     * m_class_width / 2.0;
                    double dAvrg      = Abs( from.iCno + to.iCno + 1 ) * m_class_width / 2.0 + m_range_min;
                    double dCounts    = 0.5;
                    
                    YieldDamageOverTurningPoints( from.nIdx_TP, to.nIdx_TP, dAmplitude, dAvrg, dCounts );
                    
                    from = to;
                } /* end for */
            } /* end if */
            break;
            
        case RESIDUUM_FULLCYCLES:
            if( Residuum.size() >= 2 ) 
            {
                typename VectorResiduum::const_iterator It;
                t_res from, to;
                
                from = Residuum[0];
                for ( It = Residuum.begin() + 1; It != Residuum.end(); It++ ) 
                {
                    to = *It;
                    AddCycles( from.iCno, to.iCno, 1 );
                    
                    double dAmplitude = Abs( from.iCno - to.iCno )     * m_class_width / 2.0;
                    double dAvrg      = Abs( from.iCno + to.iCno + 1 ) * m_class_width / 2.0 + m_range_min;
                    
                    YieldDamageOverTurningPoints( from.nIdx_TP, to.nIdx_TP, dAmplitude, dAvrg );
                    
                    from = to;
                } /* end for */
            } /* end if */
            break;
            
        case RESIDUUM_CLORMANN_SEEGER:
            // Zaehlt zusaetzlich Hysteresen aus Memory 1 Verhalten mit wechselndem Vorzeichen
            if( Residuum.size() >= 3 ) 
            {
                VectorResiduum res, stack;
                t_res first;
                size_t i, tail;

                // Die tatsaechlichen Werte in first sind irrelevant, da dieser nachfolgend niemals verrechnet wird
                first.iCno      = ClassNo( static_cast<T>( 0.0 ) );
                first.ClassMean = ( static_cast<T>( (double)first.iCno ) + 0.5 ) * m_class_width + m_range_min;

                res.push_back( first );
                res.insert( res.end(), Residuum.begin(), Residuum.end() );
                
                for ( i = 0; i < res.size(); i++ ) 
                {
                    bool do_loop;
                    stack.push_back( res[ i ] );
                    tail = stack.size();

                    do 
                    {
                        do_loop = tail >= 4;
                        if( do_loop ) 
                        {
                            double y1 = stack[ tail - 4 ].ClassMean;
                            double y2 = stack[ tail - 3 ].ClassMean;
                            double y3 = stack[ tail - 2 ].ClassMean;
                            double y4 = stack[ tail - 1 ].ClassMean;
                            t_res y4s = stack[ tail - 1 ];

                            do_loop = y2 * y3 < 0 && Abs( y4 ) >= Abs( y2 ) && Abs( y2 ) >= Abs( y3 );
                            if( do_loop ) 
                            {
                                t_res from = stack[ tail - 3 ]; /* y2 */
                                t_res to   = stack[ tail - 2 ]; /* y3 */
                                AddCycles( from.iCno, to.iCno, 1 );
                    
                                //size_t idx;
                                double dAmplitude  = Abs( from.iCno - to.iCno )     * m_class_width / 2.0;
                                double dAvrg       = Abs( from.iCno + to.iCno + 1 ) * m_class_width / 2.0 + m_range_min;
                    
                                YieldDamageOverTurningPoints( from.nIdx_TP, to.nIdx_TP, dAmplitude, dAvrg );
                    
                                stack[ tail - 3 ] = y4s;
                                stack.pop_back();
                                stack.pop_back();
                                tail -= 2;
                            } /* end if */
                        } /* end if */
                    } while( do_loop );
                } /* end for */
            } /* end if */
            break;
                
        case RESIDUUM_REPEATED:
            /* "simplified version" */
            if ( Residuum.size() >= 2 ) 
            {
                VectorData stream;
                
                t_matrix_backup backup;
                int* residuum_matrix;

                /* Wiederholter Durchlauf: Das Residuum wird 2x
                 * hintereinander kopiert */
                for( size_t n = 0; n < Residuum.size(); n++ ) 
                {
                    stream.push_back( Residuum[n].ClassMean );
                } /* end for */
                stream.insert( stream.end(), stream.begin(), stream.end() ); 

                m_TurningPoints[1].swap( m_TurningPoints[0] );
                m_Residuum[1].swap( m_Residuum[0] );
                // m_TurningPoints[0] und Residuum[0] wird von DoRainflow() geleert
                
                /* Klassierung durchfuehren */
                PushMatrix( &backup );
                DoRainflow( &stream[0], stream.size() );

                residuum_matrix = m_rf_matrix_no_res;
                m_rf_matrix_no_res = NULL;  // Vor dem RestoreMatrix() in Sicherheit bringen

                ASSERT( residuum_matrix != NULL );

                RestoreMatrix( &backup );
                m_TurningPoints[1].swap( m_TurningPoints[0] );
                m_Residuum[1].swap( m_Residuum[0] );

                /* Die zusaetzlichen Punkte in die 
                 * Matrix uebernehmen */
                for( int i = 0; i < m_class_count; i++ ) 
                {
                    for( int j = 0; j < m_class_count; j++ ) 
                    {
                        m_rf_matrix[ i * m_class_count + j ] += residuum_matrix[ i * m_class_count + j ];
                    } // end for
                } // end for

                free( residuum_matrix );
            } /* end if */
            break;

        case RESIDUUM_RP_DIN:
            /* Aequivalenz zu "Spannenpaar nach DIN 45667" */
            if( Residuum.size() > 2 ) 
            {
                size_t n, up_len = 0, down_len = 0;
                std::vector<t_slope> up, down;
                
                for( n = 1; n < Residuum.size(); n++ ) 
                {
                    t_res from = Residuum[n-1];
                    t_res to   = Residuum[n];
                    t_slope slope;

                    slope.dAmplitude = Abs( from.ClassMean - to.ClassMean );
                    slope.from = from;
                    slope.to = to;
                    
                    if( from.iCno > to.iCno ) 
                    {
                        /* Absteigende Aeste */
                        down.push_back( slope );
                    } else {
                        /* Aufsteigende Aeste */
                        up.push_back( slope );
                    } /* end if */
                } /* end for */
                
                /* Aufsteigende Aeste nach Amplitude sortieren */
                std::sort( up.begin(), up.end() );

                /* Absteigende Aeste nach Amplitude sortieren */
                std::sort( down.begin(), down.end() ); 
                
                /* amplitude stripping */
                /* Ref: LMS TecWare Manual, Volume 2 - Analysis, Section 4.2.4, Page 315 */
                for( n = 0; n < MIN( up.size(), down.size() ); n++ ) 
                {
                    t_slope temp = ( up[n].dAmplitude > down[n].dAmplitude ) ? down[n] : up[n];
                    AddCycles( temp.from.iCno, temp.to.iCno, 1 );

                    double dAmplitude = Abs( temp.from.iCno - temp.to.iCno ) * m_class_width / 2.0;
                    double dAvrg      = Abs( temp.from.iCno + temp.to.iCno + 1 ) * m_class_width / 2.0 + m_range_min;

                    YieldDamageOverTurningPoints( temp.from.nIdx_TP, temp.to.nIdx_TP, dAmplitude, DBL_NAN );
                } /* end for */
            } /* end if */
            break;

        default:
            ASSERT( false ); 
            break;

    } /* end switch */
} // end of DoRainflow


template<class T> 
double CRainflowT<T>::CalcDamage( double dAmplitude, double dAvrg, double dCounts ) const
{
    double dDamage;
    double dLW;
    AmplTransform AmplTransform;


    switch( m_AmplTrans_mode )
    {
        case 1:
            dAmplitude = AmplTransform.Calc( dAmplitude, dAvrg, m_AmplTrans_dM, m_AmplTrans_dR, true /* bZielIstR */ );
            break;

        case 2:
            dAmplitude = AmplTransform.Calc( dAmplitude, dAvrg, m_AmplTrans_dM, m_AmplTrans_dAvrg, false /* bZielIstR */ );   
    }

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
