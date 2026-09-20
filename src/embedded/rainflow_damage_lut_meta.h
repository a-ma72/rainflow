/**
 * @file rainflow_damage_lut_meta.h
 * @brief Single source of truth for the Woehler exponent k used to
 *        generate the currently linked RF_DAMAGE_LUT (rainflow_config.c).
 *
 * Background: a runtime correction of the form
 *
 *   damage_corrected = damage_double * (RF_CLASS_WIDTH/2/SD)^k / ND
 *
 * and the S-N curve already baked into RF_DAMAGE_LUT class-wise as
 * damage(amplitude) = (amplitude/SD)^k must use the same exponent k —
 * otherwise two different S-N curves are mixed and the combined damage
 * sum is inconsistent (looks plausible, but is wrong).
 *
 * When RF_DAMAGE_LUT is regenerated (generate_rainflow_config.py), this
 * value MUST be kept in sync — ideally the LUT generator writes it
 * out directly instead of updating it here by hand.
 *
 * PLACEHOLDER STATUS: RF_DAMAGE_LUT in rainflow_config.c is currently
 * the identity (LUT[i] = i), NOT actually computed from
 * (amplitude/SD)^k (see TODO there). RF_DAMAGE_LUT_K is still provided
 * so callers can use the same exponent and the later real generator
 * has somewhere to write it.
 */
#ifndef RAINFLOW_DAMAGE_LUT_META_H
#define RAINFLOW_DAMAGE_LUT_META_H

#define RF_DAMAGE_LUT_K 4.0

#endif /* RAINFLOW_DAMAGE_LUT_META_H */
