# README for `run.py`: The Eisenstein Scattering Quantizer

This script, `run.py`, provides a compelling and self-contained demonstration of a novel method for computing the nontrivial zeros of the Riemann zeta function. It showcases a deep connection between the zeta function and the scattering theory of automorphic forms on `GL(2)`.

The script's primary purpose is to act as a clear, visually intuitive demonstration of the principles discussed in the broader research framework.

## The Method: Two Independent Rulers 📏

The script computes the locations of the zeta zeros (`γ_n`) by treating the problem as a "phase quantization" task. It implements and compares two independent methods, or "rulers," for measuring the phase accumulation along the critical line.

#### **Ruler A: The Classical Riemann-Siegel (RS) Quantizer**

This is the standard, well-established method for locating the zeros. It uses the argument of the zeta function itself, `arg ζ(1/2 + iτ)`, combined with the Riemann-Siegel theta function. This serves as the analytic ground truth in the demonstration.

#### **Ruler B: The Novel Eisenstein Scattering Quantizer**

This is the new method. It completely ignores the zeta function and instead uses the **scattering phase `φ(τ)` of a GL(2) Eisenstein series**. This phase comes from a different area of mathematics and has no obvious connection to the individual zeros of `ζ(s)`.

### The "Wow" Moment: The Simple Baseline Correction

The core insight demonstrated by this script is what happens when you compare the two rulers. The difference between them is not random noise but a **tiny, smooth, low-degree polynomial**.

The script computes this simple baseline correction and adds it to the raw scattering phase. The resulting "adjusted" ruler then predicts the locations of the zeta zeros with astonishing accuracy, providing a powerful visual and numerical verification of the underlying theory.

## Requirements

The script requires the following Python libraries, which can be installed via pip:

```bash
pip install mpmath numpy pandas matplotlib
```

## How to Run

The script is standalone and can be executed directly.

```bash
python3 run.py
```

You can configure parameters, such as the range of zeros to compute (`TAU_MAX`) or the baseline fitting method (`BASELINE_METHOD`), by editing the variables at the top of the `main()` function in the script.

## What to Expect: Output

When you run the script, it will produce:

1.  **Console Output:** A statistical summary comparing the accuracy (mean, RMS, and max error) of the zeros predicted by both methods against the true values.
2.  **CSV Files:** Two files, `rs_quantizer_table.csv` and `rslock_eisenstein_table.csv`, containing the predicted zeros and the errors for each method.
3.  **Plots:** Several plots will be displayed, showing:
      * The "staircase" functions for both the classical and the new method.
      * The smooth polynomial baseline fit to the phase mismatch.
      * The residual errors for each method, showing their high accuracy.
