#!/usr/bin/env python3
"""
PK Parameter Calculator for PK-Sim simulation results.

Usage:
  python calculate_pk.py --input results.csv --output pk_params.json
  python calculate_pk.py --input results.csv --dose 100 --dose-unit mg

Input CSV format:
  Time [h], Concentration [mg/L]
  0, 0
  0.5, 5.2
  ...

Output: JSON with Cmax, Tmax, AUC0_t, AUC0_inf, t_half, CL, Vd, Vss, MRT
"""

import argparse
import csv
import json
import math
import sys
from typing import List, Tuple, Optional


def read_csv(filepath: str, time_col: int = 0, conc_col: int = 1) -> List[Tuple[float, float]]:
    """Read time-concentration data from CSV."""
    data = []
    with open(filepath, "r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) <= max(time_col, conc_col):
                continue
            try:
                t = float(row[time_col].strip())
                c = float(row[conc_col].strip())
                data.append((t, c))
            except ValueError:
                continue  # skip header or non-numeric rows
    return sorted(data, key=lambda x: x[0])


def trapezoidal_auc(data: List[Tuple[float, float]]) -> float:
    """Calculate AUC0-t using linear trapezoidal rule."""
    auc = 0.0
    for i in range(1, len(data)):
        dt = data[i][0] - data[i - 1][0]
        if dt <= 0:
            continue
        auc += 0.5 * (data[i][1] + data[i - 1][1]) * dt
    return auc


def log_trapezoidal_auc(data: List[Tuple[float, float]]) -> float:
    """Calculate AUC0-t using log-linear trapezoidal rule (better for elimination phase)."""
    auc = 0.0
    for i in range(1, len(data)):
        dt = data[i][0] - data[i - 1][0]
        if dt <= 0:
            continue
        c1, c2 = data[i - 1][1], data[i][1]
        if c1 > 0 and c2 > 0:
            # Log-linear interpolation
            auc += dt * (c2 - c1) / math.log(c2 / c1)
        else:
            # Fall back to linear
            auc += 0.5 * (c1 + c2) * dt
    return auc


def estimate_lambdaz(data: List[Tuple[float, float]], min_points: int = 3) -> Tuple[float, float, List[int]]:
    """Estimate terminal elimination rate constant (lambda_z) from log-linear regression.
    
    Automatically selects the best terminal phase (maximum adjusted R^2).
    Returns: (lambdaz, r_squared, indices_used)
    """
    best_r2 = -float("inf")
    best_lambdaz = 0.0
    best_indices = []
    
    n = len(data)
    # Start from the last point and work backwards
    for start_idx in range(n - min_points):
        subset = data[start_idx:]
        times = [d[0] for d in subset]
        lnc = [math.log(d[1]) for d in subset if d[1] > 0]
        t_sub = [d[0] for d, c in zip(subset, [d[1] for d in subset]) if c > 0]
        
        if len(t_sub) < min_points:
            continue
        
        # Linear regression: ln(C) = ln(C0) - lambdaz * t
        n_pts = len(t_sub)
        sum_t = sum(t_sub)
        sum_lnc = sum(lnc)
        sum_t2 = sum(t * t for t in t_sub)
        sum_tlnc = sum(t * lc for t, lc in zip(t_sub, lnc))
        
        denom = n_pts * sum_t2 - sum_t ** 2
        if denom == 0:
            continue
        
        slope = (n_pts * sum_tlnc - sum_t * sum_lnc) / denom
        intercept = (sum_lnc - slope * sum_t) / n_pts
        
        # R^2
        ss_res = sum((lc - (intercept + slope * t)) ** 2 for t, lc in zip(t_sub, lnc))
        ss_tot = sum((lc - sum_lnc / n_pts) ** 2 for lc in lnc)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        
        # Adjusted R^2
        adj_r2 = 1 - (1 - r2) * (n_pts - 1) / (n_pts - 2) if n_pts > 2 else r2
        
        if adj_r2 > best_r2 and slope < 0:
            best_r2 = adj_r2
            best_lambdaz = -slope
            best_indices = list(range(start_idx, n))
    
    return best_lambdaz, best_r2, best_indices


def calculate_pk(data: List[Tuple[float, float]], dose: Optional[float] = None) -> dict:
    """Calculate non-compartmental PK parameters."""
    if not data or len(data) < 2:
        return {"error": "Insufficient data points"}
    
    concentrations = [d[1] for d in data]
    times = [d[0] for d in data]
    
    # Cmax and Tmax
    cmax = max(concentrations)
    tmax = times[concentrations.index(cmax)]
    
    # AUC0-t (linear trapezoidal)
    auc0_t = trapezoidal_auc(data)
    
    # Lambda_z estimation
    lambdaz, r2_lambdaz, terminal_indices = estimate_lambdaz(data)
    
    # AUC0-inf
    ct = data[-1][1]
    auc0_inf = None
    if lambdaz > 0 and ct > 0:
        auc0_inf = auc0_t + ct / lambdaz
    
    # Half-life
    t_half = math.log(2) / lambdaz if lambdaz > 0 else None
    
    # Clearance and Volume
    cl = None
    vd = None
    vss = None
    mrt = None
    
    if dose is not None and auc0_inf is not None and auc0_inf > 0:
        cl = dose / auc0_inf
        vd = cl / lambdaz if lambdaz > 0 else None
        
        # AUMC calculation (first moment)
        aumc = 0.0
        for i in range(1, len(data)):
            dt = data[i][0] - data[i - 1][0]
            if dt <= 0:
                continue
            # t*C for AUMC using trapezoidal
            tc1 = data[i - 1][0] * data[i - 1][1]
            tc2 = data[i][0] * data[i][1]
            aumc += 0.5 * (tc1 + tc2) * dt
        
        if auc0_inf > 0 and lambdaz > 0 and ct > 0:
            # Add tail for AUMC
            aumc_inf = aumc + ct * (data[-1][0] / lambdaz + 1 / lambdaz ** 2)
            mrt = aumc_inf / auc0_inf
            vss = cl * mrt
    
    result = {
        "Cmax": round(cmax, 6),
        "Tmax": round(tmax, 4),
        "AUC0_t": round(auc0_t, 6),
        "AUC0_inf": round(auc0_inf, 6) if auc0_inf is not None else None,
        "lambda_z": round(lambdaz, 6) if lambdaz > 0 else None,
        "lambda_z_R2": round(r2_lambdaz, 6),
        "t_half": round(t_half, 4) if t_half is not None else None,
        "CL": round(cl, 6) if cl is not None else None,
        "Vd": round(vd, 4) if vd is not None else None,
        "Vss": round(vss, 4) if vss is not None else None,
        "MRT": round(mrt, 4) if mrt is not None else None,
        "N_points": len(data),
        "N_terminal": len(terminal_indices),
    }
    
    if dose is not None:
        result["Dose"] = dose
    
    return result


def main():
    parser = argparse.ArgumentParser(description="Non-compartmental PK parameter calculator")
    parser.add_argument("--input", required=True, help="Input CSV file (Time, Concentration)")
    parser.add_argument("--output", default=None, help="Output JSON file (default: stdout)")
    parser.add_argument("--dose", type=float, default=None, help="Dose in mg (required for CL, Vd, Vss)")
    parser.add_argument("--time-col", type=int, default=0, help="Time column index (default: 0)")
    parser.add_argument("--conc-col", type=int, default=1, help="Concentration column index (default: 1)")
    args = parser.parse_args()
    
    data = read_csv(args.input, args.time_col, args.conc_col)
    if not data:
        print("Error: No valid data found in input file", file=sys.stderr)
        sys.exit(1)
    
    pk = calculate_pk(data, args.dose)
    
    output = json.dumps(pk, indent=2)
    if args.output:
        with open(args.output, "w") as f:
            f.write(output)
        print(f"PK parameters written to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
