#Explicit Electrothermal LLP VO2 Model Reproducing Preisach Like Hysteresis for
#Memristive and Neuromorphic Devices

#Scientific Reports
#B. A. S. F. Sena 1,2,+,* and L. A. L. de Almeida1,+

#1 Federal University of ABC (UFABC), Center for Engineering, Modeling and
#Applied Social Sciences, Santo Andre, SP, 09210-580, Brazil

#2 Federal Institute of Sao Paulo (IFSP), Department of Electrical,
#São Paulo, SP, 01109-010, Brazil ˜

#* sena.bruno@ifsp.edu.br
#+ these authors contributed equally to this work


#This script reproduces
#Figure 5 from the main article
#FORC protocol and corresponding hysteretic response


#This script was originally executed in a Google Colab environment;
#running it there is recommended for faster and smoother execution.


import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os

print(f"--- O script está rodando a partir de: {os.getcwd()} ---")

# ==============================================================================
# 1. Temperature Protocol Generator Class
# ==============================================================================

class TriangularTemperatureProtocol:
    """

    """

    def __init__(self, T_min=20, T_max=80, n=10, tau_0=10, P_k=None,
                 v_up=None, v_down=None, use_linear=False):
        self.T_min = T_min
        self.T_max = T_max
        self.n = int(n)
        self.Delta_T = (T_max - T_min) / self.n
        self.tau_0 = tau_0
        self.use_linear = use_linear

        if P_k is None: self.P_k = np.ones(self.n) * 20
        elif np.isscalar(P_k): self.P_k = np.ones(self.n) * P_k
        else: self.P_k = np.array(P_k)

        if use_linear:
            self.v_up = v_up if v_up is not None else 3.0
            self.v_down = v_down if v_down is not None else 3.0

        self.a = np.zeros(self.n + 1)
        self.a[1] = tau_0
        for k in range(1, self.n):
            if use_linear:
                duration_k = (k * self.Delta_T) / self.v_down + (k * self.Delta_T) / self.v_up
                self.a[k + 1] = self.a[k] + duration_k
            else:
                self.a[k + 1] = self.a[k] + self.P_k[k-1]

        self.T_k = T_max - np.arange(1, self.n + 1) * self.Delta_T

    def temperature(self, t):
        t = np.atleast_1d(t)
        T = np.zeros_like(t, dtype=float)

        mask_init = (t >= 0) & (t < self.tau_0)
        T[mask_init] = self.T_min + (self.T_max - self.T_min) / self.tau_0 * t[mask_init]

        for k in range(1, self.n + 1):
            a_k = self.a[k]
            if self.use_linear:
                tau_down = (k * self.Delta_T) / self.v_down
                b_k = a_k + tau_down
                tau_up = (k * self.Delta_T) / self.v_up
                c_k = b_k + tau_up
                mask_down = (t >= a_k) & (t < b_k)
                T[mask_down] = self.T_max - self.v_down * (t[mask_down] - a_k)
                mask_up = (t >= b_k) & (t < c_k)
                T[mask_up] = self.T_k[k-1] + self.v_up * (t[mask_up] - b_k)
            else:
                P_k = self.P_k[k-1]
                b_k = a_k + P_k
                mask_k = (t >= a_k) & (t < b_k)
                if np.any(mask_k):
                    amplitude = self.T_max - self.T_k[k-1]
                    T[mask_k] = self.T_max - amplitude * self.tri((t[mask_k] - a_k) / P_k)


        if t.size > 0 and abs(t[-1] - self.get_total_duration()) < 1e-9:
            T[-1] = self.T_max

        return T

    def get_reversal_points(self):
        reversal_times = [0, self.tau_0]
        reversal_temps = [self.T_min, self.T_max]
        for k in range(1, self.n + 1):
            a_k = self.a[k]
            if self.use_linear:
                tau_down = (k * self.Delta_T) / self.v_down
                b_k = a_k + tau_down
                tau_up = (k * self.Delta_T) / self.v_up
                c_k = b_k + tau_up
                reversal_times.extend([b_k, c_k])
                reversal_temps.extend([self.T_k[k-1], self.T_max])
            else:
                P_k = self.P_k[k-1]
                reversal_times.extend([a_k + P_k/2, a_k + P_k])
                reversal_temps.extend([self.T_k[k-1], self.T_max])
        return np.array(reversal_times), np.array(reversal_temps)

    def get_total_duration(self):
        rev_t, _ = self.get_reversal_points()
        return rev_t[-1]

    def generate_sampled_data(self, sample_time_s=0.1):
        if sample_time_s <= 0: raise ValueError("O tempo de amostragem deve ser positivo.")
        t_total = self.get_total_duration()
        num_points = round(t_total / sample_time_s) + 1
        time_array = np.linspace(0, t_total, int(num_points))
        temperature_array = self.temperature(time_array)
        return time_array, temperature_array

# ==============================================================================
# 2. Hysteresis Model
# ==============================================================================
def ProFuncVO2(x, gama):
    return 0.5*(1 - np.sin(gama*x))*(1 + np.tanh(np.pi**2 - 2*np.pi*x))
def anchor_major_auto(T0, T1, w, Tc, beta, gama, Tpr0=None, eps_g=1e-6):
    dT = T1 - T0
    delta = +1 if dT > 0 else -1 if dT < 0 else +1
    Tr = float(T0)
    if Tpr0 is None: Tpr0 = float(w)
    P0 = ProFuncVO2(0.0, gama)
    arg0 = delta*w/2 + Tc - (Tr + Tpr0*P0)
    gr = 0.5 + 0.5*np.tanh(beta*arg0)
    gr = float(np.clip(gr, eps_g, 1.0 - eps_g))
    Tpr = float(Tpr0)
    return Tr, gr, Tpr, delta
def update_hysteresis_given_T(T, w, Tc, beta, gama, Tpr0=None, eps_dT=1e-9):
    N = len(T)
    g = np.zeros(N, dtype=float); delta_vec = np.zeros(N, dtype=int)
    Tr_vec = np.zeros(N, dtype=float); Tpr_vec = np.zeros(N, dtype=float)
    Tr, gr, Tpr, delta = anchor_major_auto(T[0], T[1], w, Tc, beta, gama, Tpr0=Tpr0)
    g[0] = gr; delta_vec[0] = delta; Tr_vec[0] = Tr; Tpr_vec[0] = Tpr
    for n in range(1, N):
        dT = T[n] - T[n-1]
        delta_new = +1 if dT >  eps_dT else (-1 if dT < -eps_dT else delta)
        if delta_new != delta:
            Tr = T[n-1]; gr = float(np.clip(g[n-1], 1e-6, 1.0 - 1e-6)); delta = delta_new
            P0 = ProFuncVO2(0.0, gama)
            Tpr = (delta*w/2 + Tc - Tr - np.arctanh(2*gr - 1)/beta) / max(P0, 1e-9)
            if abs(Tpr) < 1e-9: Tpr = np.sign(Tpr)*1e-9 if Tpr != 0 else 1e-3
        Tpr_safe = np.sign(Tpr) * max(abs(Tpr), 1e-9)
        x   = (T[n] - Tr) / Tpr_safe
        arg = delta*w/2 + Tc - (T[n] + Tpr_safe*ProFuncVO2(x, gama))
        g[n] = 0.5 + 0.5*np.tanh(beta*arg); g[n] = float(np.clip(g[n], 0.0, 1.0))
        delta_vec[n] = delta; Tr_vec[n] = Tr; Tpr_vec[n] = Tpr
    return g, delta_vec, Tr_vec, Tpr_vec
def R_of_T(T, g, Rs=17.0, Rm=140.0):
    R = g*Rs*np.exp(2553.0/(T + 273.0)) + Rm
    return np.maximum(R, 1e-12)

# ==============================================================================
# 3. Protocol Constructor
# ==============================================================================
def build_FORC_reset_protocol(Tmin=20.0, Tmax=80.0, n_levels=10, v_cool=8.0, v_heat=8.0, Fs=100.0):
    tau0 = (Tmax - Tmin) / v_heat
    protocol = TriangularTemperatureProtocol(
        T_min=Tmin, T_max=Tmax, n=n_levels,
        tau_0=tau0, v_down=v_cool, v_up=v_heat, use_linear=True
    )
    sample_time = 1.0 / Fs
    t, T = protocol.generate_sampled_data(sample_time_s=sample_time)
    rev_t, _ = protocol.get_reversal_points()
    idx_levels = []
    for k in range(1, n_levels + 1):
        t_start_level = rev_t[2*k - 1]
        t_end_level = rev_t[2*k + 1]
        i_start = np.searchsorted(t, t_start_level, side='left')
        i_end = np.searchsorted(t, t_end_level, side='right') - 1
        i_end = max(i_start, i_end)
        idx_levels.append((i_start, i_end))
    return t, T, idx_levels, protocol.Delta_T

# ==============================================================================
# 4. Main Execution and Plotting Script
# ==============================================================================
def run_and_plot():
    w, Tc, beta, gama = 6.5, 47.6, 0.2, 0.9
    Tmin, Tmax, n_levels = 30.0, 80.0, 22

    t, T, idx_levels, dT = build_FORC_reset_protocol(Tmin=Tmin, Tmax=Tmax, n_levels=n_levels, v_cool=10.0, v_heat=10.0)

    print("Simulando o modelo de histerese...")
    g, delta_vec, Tr_vec, Tpr_vec = update_hysteresis_given_T(T, w, Tc, beta, gama, Tpr0=w)
    R = R_of_T(T, g)
    print("Simulação concluída.")


    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12.5, 3.2))

    # (a) Temperature vs Time Protocol
    ax1.plot(t, T,color='black', linewidth=1.0)
    ax1.set_xlabel('Time (s)', fontsize=9)
    ax1.set_ylabel('Temperature (°C)', fontsize=9)
    ax1.set_title(f'(a) FORC protocol — {n_levels} levels', fontsize=9)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.tick_params(labelsize=8)

    # (b) FORCs g × T
    ax2.plot(T[:idx_levels[0][0]], g[:idx_levels[0][0]],
         color='gray', linewidth=0.8, label='Reset')

    for level, (i0, i1) in enumerate(idx_levels, start=1):

    # Low levels = strong, high levels = transparent.
        alpha = 1.0 - 0.75*(level-1)/(n_levels-1)

        ax2.plot(
            T[i0:i1+1],
            g[i0:i1+1],
            color='black',
            alpha=alpha,
            linewidth=1.2
        )

    ax2.set_xlabel('Temperature (°C)', fontsize=9)
    ax2.set_ylabel('Volumetric Fraction - g', fontsize=9)
    ax2.set_title('(b) FORCs g × T', fontsize=9)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.tick_params(labelsize=8)

    # Indication of the FORC level ordering
    ax2.text(
        0.97, 0.97,
        'k = 22 (highest $T_r$)\n⋮\nk = 1 (lowest $T_r$)',
        transform=ax2.transAxes,
        fontsize=8,
        verticalalignment='top',
        horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )


    if n_levels <= 10:
        ax2.legend(fontsize=6, ncol=2, loc='best')

    # (c) R vs T (Logarithmic scale)
    ax3.plot(T[:idx_levels[0][0]], R[:idx_levels[0][0]], color='black', linewidth=0.8, label='Reset')
    for i, (i0, i1) in enumerate(idx_levels, start=1):
        ax3.plot(T[i0:i1+1], R[i0:i1+1], color='black', linewidth=0.8, label=f'N{i}' if n_levels <= 10 else None)
    ax3.set_xlabel('Temperature (°C)', fontsize=9)
    ax3.set_ylabel('Resistance (Ω)', fontsize=9)
    ax3.set_yscale('log')
    ax3.set_title(f'(c) R(T) — {n_levels} Levels', fontsize=9)
    ax3.grid(True, which="both", linestyle='--', alpha=0.6)
    ax3.tick_params(labelsize=8)

    # Ordering of FORC levels
    ax3.text(
        0.97, 0.97,
        'k = 22 (highest $T_r$)\n⋮\nk = 1 (lowest $T_r$)',
        transform=ax3.transAxes,
        fontsize=8,
        verticalalignment='top',
        horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )

    if n_levels <= 10:
        ax3.legend(fontsize=6, ncol=2, loc='best')

    #
    plt.tight_layout()


    from pathlib import Path

    outdir = Path("forc_output")
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Salvando dados de saída em: {outdir.resolve()}")

    # Salvar PNG e PDF
    fig.savefig(outdir / 'forc_figure.png', dpi=300, bbox_inches='tight')
    fig.savefig(outdir / 'forc_figure.pdf', bbox_inches='tight')
    print(f"Figura salva em: {outdir / 'forc_figure.png'}")
    print(f"Figura salva em: {outdir / 'forc_figure.pdf'}")

    # Download to the computer (Colab)
    from google.colab import files
    files.download(outdir / 'forc_figure.png')
    files.download(outdir / 'forc_figure.pdf')


    plt.show()
    np.savetxt(outdir/'trajectory_all.csv', np.c_[t, T, g, R, delta_vec], delimiter=',', header='t_s,T_C,g,R_ohm,delta', comments='')
    for i, (i0, i1) in enumerate(idx_levels, start=1):
        arr = np.c_[t[i0:i1+1], T[i0:i1+1], g[i0:i1+1], R[i0:i1+1], delta_vec[i0:i1+1]]
        np.savetxt(outdir/f'forc_level_{i:02d}.csv', arr, delimiter=',', header='t_s,T_C,g,R_ohm,delta', comments='')

if __name__ == "__main__":
    run_and_plot()
