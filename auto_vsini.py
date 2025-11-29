def estimate_continuum_savgol(wave, flux, window_angstrom=10.0):
    """
    Estima el continuo usando un filtro de Savitzky-Golay sobre una
    ventana de anchura fija en Angstrom (por defecto 10 Å).

    window_angstrom : anchura de la ventana de suavizado en Angstrom.
    """
    if len(wave) < 11:
        return flux.copy()

    dw = wave[1] - wave[0]
    if dw <= 0:
        dw = abs(dw) if dw != 0 else 1.0

    # Convertir anchura en Angstrom a número de píxeles
    window_pix = int(window_angstrom / dw)
    window_pix = max(11, window_pix)

    # Asegurar longitud impar y que quepa en el vector
    if window_pix % 2 == 0:
        window_pix += 1
    if window_pix >= len(flux):
        window_pix = len(flux) - 1 if len(flux) % 2 == 0 else len(flux)

    if window_pix < 5:
        return flux.copy()

    smooth = savgol_filter(flux, window_length=window_pix, polyorder=3)
    return smooth

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
auto_vsini.py

Modo automático para el método de Fourier:
- Detecta líneas de absorción en un espectro 1D.
- Ajusta una o dos gaussianas para cada línea.
- Descarta líneas mezcladas o de mala calidad.
- Calcula v sin i mediante fourier.Fourier para cada línea aceptada.
- Guarda los resultados en un archivo ASCII ("Fourier_auto.out" por defecto).

Este script NO requiere la interfaz gráfica (Qt).
"""

import argparse
import os
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter
from lmfit import Model
import matplotlib.pyplot as plt

import fourier  # del repositorio original


def gaussian_line(x, amp, cen, sig, slope, intercept):
    """
    Gaussiana + continuo lineal.

    amp  : amplitud (negativa para líneas en absorción)
    cen  : centro de la línea
    sig  : sigma de la gaussiana (en unidades de x, típicamente Angstrom)
    slope, intercept : término de continuo lineal.
    """
    return amp * np.exp(-(x - cen) ** 2 / (2.0 * sig ** 2)) + slope * x + intercept


def double_gaussian_line(x,
                         amp1, cen1, sig1,
                         amp2, cen2, sig2,
                         slope, intercept):
    """
    Suma de dos gaussianas + continuo lineal.
    """
    g1 = amp1 * np.exp(-(x - cen1) ** 2 / (2.0 * sig1 ** 2))
    g2 = amp2 * np.exp(-(x - cen2) ** 2 / (2.0 * sig2 ** 2))
    return g1 + g2 + slope * x + intercept


def compute_bic(result):
    """
    Calcula un BIC sencillo a partir de un resultado de lmfit.

    BIC = n * ln(RSS/n) + k * ln(n)
    donde:
      n = número de puntos de datos
      k = número de parámetros libres
      RSS = chi^2 (suma de residuos al cuadrado)
    """
    n = getattr(result, "ndata", None)
    k = getattr(result, "nvarys", None)
    if n is None:
        n = len(result.data)
    if k is None:
        k = len(result.var_names)
    rss = result.chisqr
    if n <= 0:
        n = len(result.data)
    bic = n * np.log(rss / n) + k * np.log(n)
    return bic


def classify_delta_bic(delta):
    """
    Return an English string describing the strength of ΔBIC (= BIC_single - BIC_double).
    """
    if not np.isfinite(delta):
        return "ΔBIC undefined"
    if delta < 0:
        return "single model favoured"
    if delta < 1:
        return "no clear preference"
    if delta < 3:
        return "weak evidence for double"
    if delta < 20:
        return "positive to strong evidence for double"
    return "very strong evidence for double"

# Helper: classify which model is preferred for the line itself (short label)
def classify_line_preference(delta):
    """
    Short label describing which model is preferred for the line itself.
    """
    if not np.isfinite(delta):
        return "uncertain preference"
    if delta < 0:
        return "single-line model preferred"
    if delta > 0:
        return "double-line model preferred"
    return "no clear preference"

# Helper: map internal rejection codes to human-readable English labels
def describe_rejection_reason(reason):
    """
    Map internal rejection codes to human-readable English labels.
    """
    mapping = {
        "emision": "emission-dominated profile",
        "sigma_pequena": "unresolved (very narrow) line",
        "sigma_grande": "too broad (likely continuum trend)",
        "vsini_alto": "unphysical v sin i (> 300 km/s)",
        "vsini_nan": "could not measure v sin i",
        "profundidad": "insufficient line depth",
        "BIC": "strong evidence for double component (blend)",
        "rechazo": "rejected by quality cuts",
    }
    return mapping.get(reason, reason)

def fit_single_and_double(xwin, ywin, cen0, sig0):
    """
    Ajusta una gaussiana simple y una doble gaussiana a la ventana (xwin, ywin).

    Devuelve:
      result_single, result_double, bic_single, bic_double
    """
    # --- Ajuste gaussiana simple ---
    m1 = Model(gaussian_line)
    pars1 = m1.make_params(
        amp=float(np.min(ywin) - np.median(ywin)),  # negativo
        cen=float(cen0),
        sig=float(sig0),
        slope=0.0,
        intercept=float(np.median(ywin)),
    )
    # Limitar sigma a valores razonables
    pars1["sig"].min = sig0 / 5.0
    pars1["sig"].max = sig0 * 5.0
    # Enforce absorption for single Gaussian: amplitude must be negative or zero
    pars1["amp"].max = 0.0
    result_single = m1.fit(ywin, pars1, x=xwin)

    # --- Ajuste doble gaussiana ---
    m2 = Model(double_gaussian_line)
    pars2 = m2.make_params(
        amp1=float((np.min(ywin) - np.median(ywin)) / 2.0),
        cen1=float(cen0 - sig0 / 2.0),
        sig1=float(sig0),
        amp2=float((np.min(ywin) - np.median(ywin)) / 2.0),
        cen2=float(cen0 + sig0 / 2.0),
        sig2=float(sig0),
        slope=0.0,
        intercept=float(np.median(ywin)),
    )
    pars2["sig1"].min = sig0 / 5.0
    pars2["sig1"].max = sig0 * 5.0
    pars2["sig2"].min = sig0 / 5.0
    pars2["sig2"].max = sig0 * 5.0
    # Enforce both components to be in absorption
    pars2["amp1"].max = 0.0
    pars2["amp2"].max = 0.0
    # Constrain centers of both components to stay within the local window
    pars2["cen1"].min = float(xwin.min())
    pars2["cen1"].max = float(xwin.max())
    pars2["cen2"].min = float(xwin.min())
    pars2["cen2"].max = float(xwin.max())

    result_double = m2.fit(ywin, pars2, x=xwin)

    bic_single = compute_bic(result_single)
    bic_double = compute_bic(result_double)

    return result_single, result_double, bic_single, bic_double


def is_in_telluric(region_list, lam):
    """
    Devuelve True si lam cae dentro de alguna banda tellúrica
    de region_list = [(lam_min, lam_max), ...].
    """
    for (lam_min, lam_max) in region_list:
        if lam_min <= lam <= lam_max:
            return True
    return False


def detect_lines(wave, flux, distance_pix=10, snr_sigma=3.0):
    """
    Detecta líneas en absorción como mínimos en el flujo, usando una
    estimación suave del continuo para ser sensible a líneas sobre un
    gradiente fuerte del espectro.

    distance_pix : número mínimo de píxeles entre líneas.
    snr_sigma    : umbral en sigma para considerar una detección.
    """
    # Estimar un continuo suave con Savitzky-Golay sobre 10 Å
    smooth = estimate_continuum_savgol(wave, flux, window_angstrom=10.0)

    # Residuo respecto al continuo: líneas en absorción → picos positivos
    resid = smooth - flux

    # Estimación de ruido robusto sobre el residuo
    noise = np.median(np.abs(resid - np.median(resid))) + 1e-6
    height_min = snr_sigma * noise

    peaks, props = find_peaks(
        resid,
        height=height_min,
        distance=distance_pix,
    )
    return peaks, props


def measure_vsini_for_line(wave, flux, center, width, epsilon, resolution,
                           n_boot=1001, max_retries=5):
    """
    Calcula vsini llamando fourier.Fourier muchas veces (bootstrap externo).

    center   : centro de la línea (Angstrom)
    width    : ancho total considerado (Angstrom)
    epsilon  : coeficiente de oscurecimiento al borde
    resolution : R = lambda / Delta_lambda instrumental
    """
    dlam = center / (2.0 * resolution)  # mismo criterio que en la GUI
    vsini_values = []

    for _ in range(n_boot):
        # Reintenta si fourier.Fourier lanza alguna excepción
        for _ in range(max_retries):
            try:
                v = fourier.Fourier(wave, flux, center, width, dlam, epsilon)
                if np.isfinite(v):
                    vsini_values.append(v)
                break
            except Exception:
                # Permite reintentos en caso de pequeñas inestabilidades numéricas
                continue

    if len(vsini_values) == 0:
        return np.nan, np.nan

    vsini_arr = np.asarray(vsini_values, dtype=float)
    med = np.nanmedian(vsini_arr)
    std = np.nanstd(vsini_arr)
    return med, std


def process_spectrum(
    filename,
    epsilon=0.6,
    resolution=115000.0,
    n_boot=1001,
    window_angstrom=3.5,
    bic_delta_blend=10.0,
    detect_sigma=1.5,
    min_depth_sigma=1.5,
    telluric_bands=None,
    show_plots=False,
    save_plot_path=None,
    save_grid_path=None,
):
    """
    Procesa un espectro (CSV con col1=lambda, col2=flux) y devuelve
    una lista de mediciones (una por línea aceptada).

    bic_delta_blend: si BIC_double < BIC_single - bic_delta_blend,
                     consideramos la línea "mezclada" y la descartamos.
    detect_sigma   : umbral (en sigma) para la detección inicial de líneas
                     en el residuo (más bajo → más candidatos).
    min_depth_sigma: profundidad mínima en unidades de sigma del ruido
                     local para aceptar la línea.
    """
    if telluric_bands is None:
        telluric_bands = []

    # --- Cargar espectro ---
    df = pd.read_csv(filename)
    wave = df["col1"].values.astype(float)
    flux = df["col2"].values.astype(float)

    fig = None
    ax = None
    if show_plots or save_plot_path is not None:
        fig, ax = plt.subplots()
        ax.plot(wave, flux, label="Full spectrum")
        # Estimate and show continuum Savitzky-Golay (10 Å)
        smooth_cont = estimate_continuum_savgol(wave, flux, window_angstrom=10.0)
        ax.plot(wave, smooth_cont, label="Continuum (Savitzky-Golay 10 Å)")

    # --- Detección de líneas ---
    # Distancia mínima entre líneas (más permisiva que antes)
    if len(wave) > 1:
        dw = wave[1] - wave[0]
        if dw <= 0:
            dw = abs(dw) if dw != 0 else 1.0
    else:
        dw = 1.0
    distance_pix = max(3, int(window_angstrom / dw / 4.0))
    peaks, props = detect_lines(
        wave=wave,
        flux=flux,
        distance_pix=distance_pix,
        snr_sigma=detect_sigma,
    )

    results = []
    rejected_centers = []
    accepted_fits = []
    rejected_fits = []
    half_window_pix = max(5, int(window_angstrom / dw / 2.0))

    for idx in peaks:
        # Evita bordes
        if idx < half_window_pix or idx > len(wave) - half_window_pix:
            cen0 = wave[idx]
            rejected_centers.append(cen0)
            continue

        # Ventana local alrededor del pico
        i1 = idx - half_window_pix
        i2 = idx + half_window_pix + 1
        xwin = wave[i1:i2]
        ywin = flux[i1:i2]

        cen0 = wave[idx]
        # Estimación gruesa de sigma a partir del tamaño de la ventana
        sig0 = window_angstrom / 6.0  # porque usaremos width ~ 6 sigma

        # --- Fast morphological check to reject obvious non-lines before fitting ---
        n_edge = max(3, len(ywin) // 10)
        edges = np.concatenate([ywin[:n_edge], ywin[-n_edge:]])
        cont_edge = np.median(edges)
        depth_obs = cont_edge - np.min(ywin)  # positive for absorption
        noise_edge = np.median(np.abs(edges - cont_edge)) + 1e-6
        if depth_obs < min_depth_sigma * noise_edge:
            rejected_centers.append(cen0)
            continue
        emis_ratio = (np.max(ywin) - cont_edge) / max(depth_obs, 1e-6)
        if emis_ratio > 2.0:
            rejected_centers.append(cen0)
            continue

        # Descartar si la línea está en una banda tellúrica
        if is_in_telluric(telluric_bands, cen0):
            rejected_centers.append(cen0)
            continue

        try:
            res_single, res_double, bic_s, bic_d = fit_single_and_double(
                xwin, ywin, cen0, sig0
            )
        except Exception:
            # Problemas en el ajuste → descartar
            rejected_centers.append(cen0)
            continue

        # Parámetros del ajuste simple
        best_amp = res_single.best_values["amp"]
        best_sig = res_single.best_values["sig"]
        best_cen = res_single.best_values["cen"]

        # Parámetros del ajuste doble (para diagnóstico de blends)
        amp1_d = res_double.best_values.get("amp1", np.nan)
        amp2_d = res_double.best_values.get("amp2", np.nan)
        cen1_d = res_double.best_values.get("cen1", best_cen)
        cen2_d = res_double.best_values.get("cen2", best_cen)
        sig1_d = res_double.best_values.get("sig1", best_sig)
        sig2_d = res_double.best_values.get("sig2", best_sig)
        slope2_d = res_double.best_values.get("slope", 0.0)
        intercept2_d = res_double.best_values.get("intercept", float(np.median(ywin)))

        # Rechazar explícitamente líneas en emisión (amp > 0)
        # Sólo queremos líneas en absorción para medir v sin i.
        if best_amp > 0:
            rejected_centers.append(cen0)
            if show_plots or save_grid_path is not None:
                rejected_fits.append(
                    {
                        "xwin": xwin.copy(),
                        "ywin": ywin.copy(),
                        "best_cen": best_cen,
                        "best_amp": best_amp,
                        "best_sig": best_sig,
                        "slope": res_single.best_values.get("slope", 0.0),
                        "intercept": res_single.best_values.get(
                            "intercept", float(np.median(ywin))
                        ),
                        "bic_s": bic_s,
                        "bic_d": bic_d,
                        "reason": "emision",
                        "amp1": amp1_d,
                        "amp2": amp2_d,
                        "cen1": cen1_d,
                        "cen2": cen2_d,
                        "sig1": sig1_d,
                        "sig2": sig2_d,
                        "slope2": slope2_d,
                        "intercept2": intercept2_d,
                    }
                )
            continue

        # Rechazar líneas demasiado angostas: al menos ~3 píxeles de sigma
        min_sigma = 3.0 * dw
        if best_sig < min_sigma:
            rejected_centers.append(cen0)
            if show_plots or save_grid_path is not None:
                rejected_fits.append(
                    {
                        "xwin": xwin.copy(),
                        "ywin": ywin.copy(),
                        "best_cen": best_cen,
                        "best_amp": best_amp,
                        "best_sig": best_sig,
                        "slope": res_single.best_values.get("slope", 0.0),
                        "intercept": res_single.best_values.get("intercept", float(np.median(ywin))),
                        "bic_s": bic_s,
                        "bic_d": bic_d,
                        "reason": "sigma_pequena",
                        "amp1": amp1_d,
                        "amp2": amp2_d,
                        "cen1": cen1_d,
                        "cen2": cen2_d,
                        "sig1": sig1_d,
                        "sig2": sig2_d,
                        "slope2": slope2_d,
                        "intercept2": intercept2_d,
                    }
                )
            continue

        # New: reject lines that are too broad (likely continuum, not a line)
        max_sigma = 0.3 * (xwin.max() - xwin.min())
        if best_sig > max_sigma:
            rejected_centers.append(cen0)
            if show_plots or save_grid_path is not None:
                rejected_fits.append(
                    {
                        "xwin": xwin.copy(),
                        "ywin": ywin.copy(),
                        "best_cen": best_cen,
                        "best_amp": best_amp,
                        "best_sig": best_sig,
                        "slope": res_single.best_values.get("slope", 0.0),
                        "intercept": res_single.best_values.get("intercept", float(np.median(ywin))),
                        "bic_s": bic_s,
                        "bic_d": bic_d,
                        "reason": "sigma_grande",
                        "amp1": amp1_d,
                        "amp2": amp2_d,
                        "cen1": cen1_d,
                        "cen2": cen2_d,
                        "sig1": sig1_d,
                        "sig2": sig2_d,
                        "slope2": slope2_d,
                        "intercept2": intercept2_d,
                    }
                )
            continue

        # Criterio "inteligente" para considerar la línea mezclada:
        # 1) El modelo doble mejora mucho el BIC (ΔBIC grande)
        # 2) Los dos centros están claramente separados
        # 3) El segundo componente tiene una amplitud comparable al primero
        #
        # Esto hace el criterio más conservador: sólo descartamos como
        # "mezcladas" las líneas con evidencia muy clara de doble componente.
        delta_bic = bic_s - bic_d
        amp1 = res_double.best_values.get("amp1", 0.0)
        amp2 = res_double.best_values.get("amp2", 0.0)
        cen1 = res_double.best_values.get("cen1", best_cen)
        cen2 = res_double.best_values.get("cen2", best_cen)
        sep = abs(cen2 - cen1)

        # Requerimos separación de al menos ~1 * sigma de la gaussiana simple
        min_sep = 1.0 * best_sig
        # Y que el segundo componente tenga al menos un 50% de la amplitud del primero
        amp_ratio = abs(amp2) / max(abs(amp1), 1e-6)

        is_real_blend = (
            delta_bic > bic_delta_blend
            and sep > min_sep
            and amp_ratio > 0.5
        )

        if is_real_blend:
            rejected_centers.append(cen0)
            if show_plots or save_grid_path is not None:
                rejected_fits.append(
                    {
                        "xwin": xwin.copy(),
                        "ywin": ywin.copy(),
                        "best_cen": best_cen,
                        "best_amp": best_amp,
                        "best_sig": best_sig,
                        "slope": res_single.best_values.get("slope", 0.0),
                        "intercept": res_single.best_values.get(
                            "intercept", float(np.median(ywin))
                        ),
                        "bic_s": bic_s,
                        "bic_d": bic_d,
                        "reason": "BIC",
                        "amp1": amp1_d,
                        "amp2": amp2_d,
                        "cen1": cen1_d,
                        "cen2": cen2_d,
                        "sig1": sig1_d,
                        "sig2": sig2_d,
                        "slope2": slope2_d,
                        "intercept2": intercept2_d,
                    }
                )
            continue

        # Profundidad aproximada en el centro
        depth = -best_amp  # debería ser >0 para absorción
        noise_local = np.median(np.abs(ywin - np.median(ywin))) + 1e-6
        snr = depth / noise_local if noise_local > 0 else np.nan
        if depth < min_depth_sigma * noise_local:
            rejected_centers.append(cen0)
            if show_plots or save_grid_path is not None:
                rejected_fits.append(
                    {
                        "xwin": xwin.copy(),
                        "ywin": ywin.copy(),
                        "best_cen": best_cen,
                        "best_amp": best_amp,
                        "best_sig": best_sig,
                        "slope": res_single.best_values.get("slope", 0.0),
                        "intercept": res_single.best_values.get("intercept", float(np.median(ywin))),
                        "bic_s": bic_s,
                        "bic_d": bic_d,
                        "reason": "profundidad",
                        "amp1": amp1_d,
                        "amp2": amp2_d,
                        "cen1": cen1_d,
                        "cen2": cen2_d,
                        "sig1": sig1_d,
                        "sig2": sig2_d,
                        "slope2": slope2_d,
                        "intercept2": intercept2_d,
                    }
                )
            continue

        # Anchura total que usaremos en Fourier (6 sigma)
        width = 6.0 * best_sig

        # Medida de vsini por Fourier
        vsini_med, vsini_std = measure_vsini_for_line(
            wave=wave,
            flux=flux,
            center=best_cen,
            width=width,
            epsilon=epsilon,
            resolution=resolution,
            n_boot=n_boot,
        )

        if not np.isfinite(vsini_med):
            rejected_centers.append(cen0)
            if show_plots or save_grid_path is not None:
                rejected_fits.append(
                    {
                        "xwin": xwin.copy(),
                        "ywin": ywin.copy(),
                        "best_cen": best_cen,
                        "best_amp": best_amp,
                        "best_sig": best_sig,
                        "slope": res_single.best_values.get("slope", 0.0),
                        "intercept": res_single.best_values.get("intercept", float(np.median(ywin))),
                        "bic_s": bic_s,
                        "bic_d": bic_d,
                        "reason": "vsini_nan",
                        "amp1": amp1_d,
                        "amp2": amp2_d,
                        "cen1": cen1_d,
                        "cen2": cen2_d,
                        "sig1": sig1_d,
                        "sig2": sig2_d,
                        "slope2": slope2_d,
                        "intercept2": intercept2_d,
                    }
                )
            continue

        # Rechazar valores de v sin i no físicos o sospechosos (>300 km/s)
        if np.isfinite(vsini_med) and vsini_med > 300.0:
            rejected_centers.append(cen0)
            if show_plots or save_grid_path is not None:
                rejected_fits.append(
                    {
                        "xwin": xwin.copy(),
                        "ywin": ywin.copy(),
                        "best_cen": best_cen,
                        "best_amp": best_amp,
                        "best_sig": best_sig,
                        "slope": res_single.best_values.get("slope", 0.0),
                        "intercept": res_single.best_values.get("intercept", float(np.median(ywin))),
                        "bic_s": bic_s,
                        "bic_d": bic_d,
                        "reason": "vsini_alto",
                        "amp1": amp1_d,
                        "amp2": amp2_d,
                        "cen1": cen1_d,
                        "cen2": cen2_d,
                        "sig1": sig1_d,
                        "sig2": sig2_d,
                        "slope2": slope2_d,
                        "intercept2": intercept2_d,
                    }
                )
            continue

        if show_plots and ax is not None:
            xx = np.linspace(xwin.min(), xwin.max(), 500)
            model_single = res_single.eval(x=xx)
            # Dibujar el ajuste gaussiano sobre el espectro completo
            ax.plot(xx, model_single, linewidth=1.0, alpha=0.9, color="orange")
            ax.axvline(best_cen, linestyle="--", alpha=0.7, color="orange")

        if show_plots or save_grid_path is not None:
            accepted_fits.append(
                {
                    "xwin": xwin.copy(),
                    "ywin": ywin.copy(),
                    "best_cen": best_cen,
                    "vsini_med": vsini_med,
                    "vsini_std": vsini_std,
                    "best_amp": best_amp,
                    "best_sig": best_sig,
                    "slope": res_single.best_values.get("slope", 0.0),
                    "intercept": res_single.best_values.get("intercept", float(np.median(ywin))),
                    "bic_s": bic_s,
                    "bic_d": bic_d,
                    "snr": snr,
                    "redchi_single": res_single.redchi,
                    "delta_bic": delta_bic,
                }
            )

        results.append(
            {
                "file": os.path.basename(filename),
                "lambda": best_cen,
                "width": width,
                "epsilon": epsilon,
                "R": resolution,
                "vsini": vsini_med,
                "vsini_err": vsini_std,
                "snr": snr,
                "redchi_single": res_single.redchi,
                "delta_bic": delta_bic,
            }
        )

    if (show_plots or save_plot_path is not None) and ax is not None:
        # (Opcional) líneas verticales en los centros rechazados.
        # Desactivado por el momento para evitar saturar la figura.
        # for cen_rej in rejected_centers:
        #     ax.axvline(cen_rej, linestyle="--", color="red", alpha=0.5)

        # Dibujar también las gaussianas de líneas rechazadas (cuando hay ajuste),
        # pero omitiendo las muy angostas, las de v sin i > 300 km/s, y las de emisión.
        for fit in rejected_fits:
            if "best_amp" not in fit:
                continue
            reason = fit.get("reason", "rechazo")
            if reason in ("sigma_pequena", "vsini_alto", "emision"):
                continue

            xw = fit["xwin"]
            cen = fit["best_cen"]
            amp = fit["best_amp"]
            sig = fit["best_sig"]
            slope = fit.get("slope", 0.0)
            intercept = fit.get("intercept", float(np.median(xw)))
            xx = np.linspace(xw.min(), xw.max(), 400)

            if reason == "BIC" and all(
                k in fit for k in ("amp1", "amp2", "cen1", "cen2", "sig1", "sig2")
            ):
                # Para líneas marcadas como mezcladas, dibujar el modelo doble completo
                amp1 = fit["amp1"]
                amp2 = fit["amp2"]
                cen1 = fit["cen1"]
                cen2 = fit["cen2"]
                sig1 = fit["sig1"]
                sig2 = fit["sig2"]
                slope2 = fit.get("slope2", slope)
                intercept2 = fit.get("intercept2", intercept)

                model_double = double_gaussian_line(
                    xx,
                    amp1=amp1,
                    cen1=cen1,
                    sig1=sig1,
                    amp2=amp2,
                    cen2=cen2,
                    sig2=sig2,
                    slope=slope2,
                    intercept=intercept2,
                )
                ax.plot(xx, model_double, linewidth=1.0, alpha=0.9, color="red", label=None)
            else:
                # Para otros rechazos, dibujar la gaussiana simple ajustada
                model = gaussian_line(xx, amp=amp, cen=cen, sig=sig, slope=slope, intercept=intercept)
                ax.plot(xx, model, linewidth=1.0, alpha=0.8, color="red", label=None)

        ax.set_xlabel("Wavelength (Å)")
        ax.set_ylabel("Flux")
        n_single_accepted = len(accepted_fits)
        ax.set_title(
            f"{os.path.basename(filename)}: full spectrum and Gaussian fits\n"
            f"Accepted single lines: {n_single_accepted}"
        )
        ax.legend()
        fig.tight_layout()
        if save_plot_path is not None:
            fig.savefig(save_plot_path, dpi=150)
        if show_plots:
            plt.show()
        else:
            plt.close(fig)

    # Figura adicional: grid de líneas bien ajustadas
    if (show_plots or save_grid_path is not None) and len(accepted_fits) > 0:
        # To avoid creating an excessively tall figure (which can crash matplotlib
        # with "Image size ... is too large"), limit the number of panels we plot.
        max_panels = 45  # e.g., 15 rows x 3 columns
        nlines_total = len(accepted_fits)
        nlines_plot = min(nlines_total, max_panels)
        ncols = 3
        nrows = int(np.ceil(nlines_plot / ncols))

        fig_grid, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(4 * ncols, 3 * nrows),
            squeeze=False,
        )
        axes = axes.ravel()

        if nlines_total > max_panels:
            print(
                f"[auto_vsini] Warning: {nlines_total} accepted lines, "
                f"only plotting the first {max_panels} in the grid figure to avoid an oversized image."
            )

        for i, fit in enumerate(accepted_fits[:nlines_plot]):
            axg = axes[i]
            xw = fit["xwin"]
            yw = fit["ywin"]
            cen = fit["best_cen"]
            vsmed = fit["vsini_med"]
            vsstd = fit["vsini_std"]
            amp = fit["best_amp"]
            sig = fit["best_sig"]
            slope = fit["slope"]
            intercept = fit["intercept"]
            bic_s = fit.get("bic_s", np.nan)
            bic_d = fit.get("bic_d", np.nan)
            delta_bic = bic_s - bic_d
            bic_label = classify_delta_bic(delta_bic)
            line_pref = classify_line_preference(delta_bic)
            snr = fit.get("snr", np.nan)

            xx = np.linspace(xw.min(), xw.max(), 400)
            model = gaussian_line(xx, amp=amp, cen=cen, sig=sig, slope=slope, intercept=intercept)

            axg.plot(xw, yw, label="Local spectrum")
            axg.plot(xx, model, label="Single-Gaussian model", color="orange")
            axg.axvline(cen, linestyle="--", color="orange", alpha=0.7)
            axg.set_xlabel("Wavelength (Å)")
            axg.set_ylabel("Flux")
            axg.set_title(f"λ≈{cen:.2f} Å", fontsize=8)
            info_text = (
                f"v sin i = {vsmed:.1f} ± {vsstd:.1f} km/s\n"
                f"BIC_single = {bic_s:.1f}, BIC_double = {bic_d:.1f}\n"
                f"ΔBIC = {delta_bic:.1f}  ({bic_label})\n"
                f"{line_pref}\n"
                f"depth S/N ≈ {snr:.1f}"
            )
            axg.text(
                0.03,
                0.97,
                info_text,
                transform=axg.transAxes,
                va="top",
                ha="left",
                fontsize=7,
            )
            axg.tick_params(labelsize=8)

        # Turn off unused axes
        for j in range(nlines_plot, len(axes)):
            axes[j].axis("off")

        fig_grid.tight_layout()
        if save_grid_path is not None:
            fig_grid.savefig(save_grid_path, dpi=150)
        if show_plots:
            plt.show()
        else:
            plt.close(fig_grid)

    # Figura adicional: grid de líneas rechazadas (con información de ΔBIC)
    if (show_plots or save_grid_path is not None) and len(rejected_fits) > 0:
        max_panels_rej = 45
        nrej_total = len(rejected_fits)
        nrej_plot = min(nrej_total, max_panels_rej)
        ncols_rej = 3
        nrows_rej = int(np.ceil(nrej_plot / ncols_rej))

        fig_rej, axes_rej = plt.subplots(
            nrows_rej,
            ncols_rej,
            figsize=(4 * ncols_rej, 3 * nrows_rej),
            squeeze=False,
        )
        axes_rej = axes_rej.ravel()

        if nrej_total > max_panels_rej:
            print(
                f"[auto_vsini] Warning: {nrej_total} rejected lines, "
                f"only plotting the first {max_panels_rej} in the rejected grid figure."
            )

        for i, fit in enumerate(rejected_fits[:nrej_plot]):
            axr = axes_rej[i]
            xw = fit["xwin"]
            yw = fit["ywin"]
            cen = fit["best_cen"]
            bic_s = fit.get("bic_s", np.nan)
            bic_d = fit.get("bic_d", np.nan)
            delta_bic = bic_s - bic_d
            bic_label = classify_delta_bic(delta_bic)
            reason = fit.get("reason", "rechazo")
            line_pref = classify_line_preference(delta_bic)
            reason_label = describe_rejection_reason(reason)

            axr.plot(xw, yw, label="Local spectrum")

            if all(k in fit for k in ("amp1", "amp2", "cen1", "cen2", "sig1", "sig2")):
                xx = np.linspace(xw.min(), xw.max(), 400)
                amp1 = fit["amp1"]
                amp2 = fit["amp2"]
                cen1 = fit["cen1"]
                cen2 = fit["cen2"]
                sig1 = fit["sig1"]
                sig2 = fit["sig2"]
                slope2 = fit.get("slope2", 0.0)
                intercept2 = fit.get("intercept2", float(np.median(xw)))
                model_double = double_gaussian_line(
                    xx,
                    amp1=amp1,
                    cen1=cen1,
                    sig1=sig1,
                    amp2=amp2,
                    cen2=cen2,
                    sig2=sig2,
                    slope=slope2,
                    intercept=intercept2,
                )
                axr.plot(xx, model_double, color="orange", linestyle="-", label="Double-Gaussian model")
                axr.axvline(cen1, linestyle=":", color="orange", alpha=0.7)
                axr.axvline(cen2, linestyle=":", color="orange", alpha=0.7)
            else:
                axr.axvline(cen, linestyle="--", color="red", alpha=0.7)

            axr.set_xlabel("Wavelength (Å)")
            axr.set_ylabel("Flux")
            axr.set_title(f"λ≈{cen:.2f} Å", fontsize=8)
            rej_text = (
                f"BIC_single = {bic_s:.1f}, BIC_double = {bic_d:.1f}\n"
                f"ΔBIC = {delta_bic:.1f} ({bic_label})\n"
                f"{line_pref}\n"
                f"Reason: {reason_label}"
            )
            axr.text(
                0.03,
                0.97,
                rej_text,
                transform=axr.transAxes,
                va="top",
                ha="left",
                fontsize=7,
            )
            axr.tick_params(labelsize=8)

        for j in range(nrej_plot, len(axes_rej)):
            axes_rej[j].axis("off")

        fig_rej.tight_layout()
        if save_grid_path is not None:
            base, ext = os.path.splitext(save_grid_path)
            fig_rej.savefig(f"{base}_rejected{ext or '.png'}", dpi=150)
        if show_plots:
            plt.show()
        else:
            plt.close(fig_rej)

    # Figura adicional: histograma de anchos de las gaussianas aceptadas vs rechazadas
    if (show_plots or save_grid_path is not None) and (len(accepted_fits) > 0 or len(rejected_fits) > 0):
        sig_accepted = [
            fit["best_sig"] for fit in accepted_fits
            if "best_sig" in fit and np.isfinite(fit["best_sig"])
        ]
        sig_rejected = [
            fit["best_sig"] for fit in rejected_fits
            if "best_sig" in fit and np.isfinite(fit["best_sig"])
        ]

        if len(sig_accepted) > 0 or len(sig_rejected) > 0:
            fig_hist, ax_hist = plt.subplots(figsize=(6, 4))
            bins = 20
            if len(sig_accepted) > 0:
                ax_hist.hist(
                    sig_accepted,
                    bins=bins,
                    alpha=0.6,
                    label="Accepted",
                    histtype="stepfilled",
                )
            if len(sig_rejected) > 0:
                ax_hist.hist(
                    sig_rejected,
                    bins=bins,
                    alpha=0.6,
                    label="Rejected",
                    histtype="step",
                )
            ax_hist.set_xlabel("Gaussian σ (wavelength units)")
            ax_hist.set_ylabel("Number of lines")
            ax_hist.legend()
            fig_hist.tight_layout()

            if save_grid_path is not None:
                base, ext = os.path.splitext(save_grid_path)
                if not ext:
                    ext = ".png"
                fig_hist.savefig(f"{base}_sigma_hist{ext}", dpi=150)

            if show_plots:
                plt.show()
            else:
                plt.close(fig_hist)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Modo automático para medir v sin i con el método de Fourier."
    )
    parser.add_argument(
        "spectra",
        nargs="+",
        help="Espectros de entrada en formato CSV (col1=lambda, col2=flux).",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=0.6,
        help="Coeficiente de oscurecimiento al borde (default: 0.6).",
    )
    parser.add_argument(
        "--R",
        type=float,
        default=115000.0,
        help="Resolución espectral R=lambda/dlambda (default: 115000).",
    )
    parser.add_argument(
        "--nboot",
        type=int,
        default=1001,
        help="Número de llamadas a fourier.Fourier por línea (default: 1001).",
    )
    parser.add_argument(
        "--window",
        type=float,
        default=2.0,
        help="Ancho de ventana en Angstrom alrededor de cada línea (default: 2.0).",
    )
    parser.add_argument(
        "--bic-delta",
        type=float,
        default=10.0,
        help=(
            "Umbral de ΔBIC para descartar líneas mezcladas: "
            "si BIC_doble < BIC_simple - ΔBIC, se descarta (default: 6.0)."
        ),
    )
    parser.add_argument(
        "--min-depth-sigma",
        type=float,
        default=3.0,
        help=(
            "Profundidad mínima de la línea en unidades de sigma "
            "del ruido local (default: 2.0)."
        ),
    )
    parser.add_argument(
        "--detect-sigma",
        type=float,
        default=3.0,
        help=(
            "Umbral en sigma para la detección inicial de líneas "
            "en el residuo (default: 1.0)."
        ),
    )
    parser.add_argument(
        "--output",
        type=str,
        default="Fourier_auto.out",
        help="Nombre del archivo de salida (default: Fourier_auto.out).",
    )
    parser.add_argument(
        "--show-plots",
        action="store_true",
        help="Mostrar gráficas del espectro y el ajuste gaussiano para cada línea aceptada.",
    )
    parser.add_argument(
        "--save-plot",
        type=str,
        default=None,
        help="Ruta de archivo para guardar la figura del espectro completo con ajustes.",
    )
    parser.add_argument(
        "--save-grid",
        type=str,
        default=None,
        help="Ruta de archivo para guardar un grid con las líneas bien ajustadas.",
    )

    args = parser.parse_args()

    all_results = []
    for spec in args.spectra:
        spec_results = process_spectrum(
            filename=spec,
            epsilon=args.epsilon,
            resolution=args.R,
            n_boot=args.nboot,
            window_angstrom=args.window,
            bic_delta_blend=args.bic_delta,
            detect_sigma=args.detect_sigma,
            min_depth_sigma=args.min_depth_sigma,
            telluric_bands=[],  # personalizar aquí si quieres enmascarar bandas
            show_plots=args.show_plots,
            save_plot_path=args.save_plot,
            save_grid_path=args.save_grid,
        )
        all_results.extend(spec_results)

    if not all_results:
        print("No se obtuvo ninguna medición de v sin i.")
        return

    # Escribir archivo ASCII tipo Fourier_auto.out
    now = datetime.utcnow()
    date_str = now.strftime("%Y-%m-%d")
    time_str = now.strftime("%H:%M:%S")

    header = "#file date time lambda width R epsilon vsini vsini_err snr redchi_single delta_bic\n"
    mode = "w"
    if os.path.exists(args.output):
        # Si ya existe, asumimos que ya tiene cabecera y sólo agregamos
        mode = "a"
        header = ""

    with open(args.output, mode) as f:
        if header:
            f.write(header)
        for row in all_results:
            f.write(
                f"{row['file']} {date_str} {time_str} "
                f"{row['lambda']:.4f} {row['width']:.4f} "
                f"{row['R']:.1f} {row['epsilon']:.3f} "
                f"{row['vsini']:.3f} {row['vsini_err']:.3f} {row['snr']:.3f} "
                f"{row['redchi_single']:.3f} {row['delta_bic']:.3f}\n"
            )

    print(f"Se escribieron {len(all_results)} líneas en '{args.output}'.")
    print("Ejemplo de filas:")
    for r in all_results[:5]:
        print(
            f"{r['file']}  lambda={r['lambda']:.4f}  "
            f"vsini={r['vsini']:.2f} ± {r['vsini_err']:.2f} km/s"
        )


if __name__ == "__main__":
    main()

r"""python auto_vsini.py spec_to_test.csv \ 
    --detect-sigma 3.0 \
    --min-depth-sigma 3.0 \
    --bic-delta 10.0 \
    --window 3.5 \
    --nboot 1001 \
    --show-plots \
    --save-plot spec_full.png \
    --save-grid spec_grid.png
"""