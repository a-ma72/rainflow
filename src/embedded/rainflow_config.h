/**
 * @file rainflow_config.h
 * @brief Application-specific configuration for the rainflow engine
 *        (rainflow.h/.c): class count, class bounds, and the damage LUT
 *        precomputed from the Woehler / S-N curve.
 *
 * PLACEHOLDER VALUES. This file is meant to be generated (e.g.
 * generate_rainflow_config.py from value range + S-N curve) — do not
 * maintain it by hand once the generator exists.
 *
 * Class count and value range were originally parameterized as powers
 * of two (via shift values) so rf_classify() in rainflow.c could
 * classify with a single right shift instead of a division — relevant
 * for targets without a (fast) hardware divider.
 */
#ifndef RAINFLOW_CONFIG_H
#define RAINFLOW_CONFIG_H

#include <stdint.h>

/* ------------------------------------------------------------------ */
/* Fixed-point configuration (stage 1: hysteresis filter)             */
/* ------------------------------------------------------------------ */

typedef int32_t RF_Value_t; /* Q28.4 fixed-point raw sample (before classification) */

/** Q-format: fractional bits in an int32_t. */
#define RF_FIXED_SHIFT 4U
#define RF_DOUBLE_TO_FIXED(value) ((RF_Value_t)(value * (1 << RF_FIXED_SHIFT)))
#define RF_FIXED_TO_DOUBLE(value) ((double)(value) / (1 << RF_FIXED_SHIFT))

/* ------------------------------------------------------------------ */
/* Classification configuration (stage 2: integer domain)             */
/* ------------------------------------------------------------------ */

/**
 * Number of class bins.
 * Also determines maximum residue size (2 * RF_NUM_CLASSES) and
 * damage-LUT size.
 */
#define RF_NUM_CLASSES   128
#define RF_CLASS_WIDTH   31.25
#define RF_CLASS_MIN    -1500.0

#define RF_CLASS_WIDTH_FIXED RF_DOUBLE_TO_FIXED(RF_CLASS_WIDTH)    /* Q28.4 class width */
#define RF_CLASS_MIN_FIXED   RF_DOUBLE_TO_FIXED(RF_CLASS_MIN)  /* Q28.4 lower bound */
#define RF_CLASS_MAX_FIXED   (RF_CLASS_MIN_FIXED + RF_CLASS_WIDTH_FIXED * RF_NUM_CLASSES)

/**
 * Damage contribution of a single closed cycle, indexed by class
 * difference |class_from - class_to| (0 .. RF_NUM_CLASSES-1).
 * Precomputed from the S-N curve: damage(amplitude) = (amplitude/SD)^k,
 * already scaled in fixed-point damage units.
 * TODO: fill with real S-N curve values (generated).
 */
extern const uint32_t RF_DAMAGE_LUT[RF_NUM_CLASSES];

#endif /* RAINFLOW_CONFIG_H */
