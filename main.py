from iapws import IAPWS97 as WSP
import streamlit as st
import pandas as pd
from bokeh.plotting import figure
import matplotlib.pyplot as plt
N_e = 545e6
p_0 = 26.5e6
t_0 = list(range(int(500), int(601), 1))
T_0 = [t+273.15 for t in t_0]
p_pp = 4.4e6
t_pp = 570
T_pp = t_pp+273.15
p_k = 3.4e3
t_pv = 280
T_pv = t_pv+273.15
delta_p_0 = 0.05*p_0
delta_p_pp = 0.08*p_pp
delta_p = 0.03*p_pp


def Calculate_G0_Gk(p0, T0, ppp, Tpp, pk, Tpv, Ne, deltap0, deltappp, deltap):
    point_0 = WSP(P=p0 * 1e-6, T=T0)
    h_0 = point_0.h

    p_0_d = p0 - deltap0
    point_p_0_d = WSP(P=p_0_d * 1e-6, h=point_0.h)

    p_1t = ppp + deltappp
    point_1t = WSP(P=p_1t * 1e-6, s=point_0.s)
    h_1t = point_1t.h

    point_pp = WSP(P=ppp * 1e-6, T=Tpp)
    s_pp = point_pp.s
    h_pp = point_pp.h

    H_01 = point_0.h - point_1t.h

    kpd_oi = 0.85
    H_i_cvd = H_01 * kpd_oi

    h_1 = point_0.h - H_i_cvd
    point_1 = WSP(P=p_1t * 1e-6, h=h_1)

    p_pp_d = ppp - deltappp
    point_pp_d = WSP(P=p_pp_d * 1e-6, h=point_pp.h)

    point_kt = WSP(P=pk * 1e-6, s=point_pp.s)

    H_02 = h_pp - point_kt.h

    kpd_oi = 0.85
    H_i_csd_cnd = H_02 * kpd_oi

    h_k = point_pp.h - H_i_csd_cnd
    point_k = WSP(P=pk * 1e-6, h=h_k)
    T_k = point_k.T
    h_k = point_k.h

    point_k_v = WSP(P=pk * 1e-6, x=0)
    s_k_v = point_k_v.s
    h_k_v = point_k_v.h
    p_pv = 1.4 * p0
    point_pv = WSP(P=p_pv * 1e-6, T=Tpv)

    ksi_pp_oo = 1 - (1 - (T_k * (s_pp - s_k_v)) / ((h_0 - h_1t) + (h_pp - h_k_v))) / (
                1 - (T_k * (s_pp - point_pv.s)) / ((h_0 - h_1t) + (h_pp - point_pv.h)))

    T_0_ = 374.2 + 273.15

    T_ = (point_pv.T - T_k) / (T_0_ - T_k)

    if T_ <= 0.636363636:
        ksi1 = -4.9655 * T_ ** 3 + 5.5547 * T_ ** 2 - 1.0808 * T_ + 0.4922
    elif 0.636363636 < T_ <= 0.736363636:
        ksi1 = -1.3855 * T_ ** 2 + 2.0774 * T_ + 0.0321
    elif 0.736363636 < T_ <= 0.827272727:
        ksi1 = -1.4152 * T_ ** 2 + 2.3287 * T_ - 0.1088
    else:
        ksi1 = 0.82

    if T_ <= 0.631818182:
        ksi2 = -1.0078 * T_ ** 3 - 0.3296 * T_ ** 2 + 1.7524 * T_ + 0.0714
    elif 0.631818182 < T_ <= 0.718181818:
        ksi2 = -2.5821 * T_ ** 2 + 3.689 * T_ - 0.4825
    elif 0.718181818 < T_ <= 0.827272727:
        ksi2 = -11.992 * T_ ** 3 + 25.812 * T_ ** 2 - 18.31 * T_ + 5.1453
    elif 0.827272727 < T_ <= 0.936363636:
        ksi2 = -11.115 * T_ ** 3 + 27.375 * T_ ** 2 - 22.574 * T_ + 7.1385
    else:
        ksi2 = 0.89

    ksi = (ksi1 + ksi2) / 2
    ksi_pp = ksi * ksi_pp_oo

    kpd_ir = ((H_i_cvd + H_i_csd_cnd) / (H_i_cvd + (h_pp - h_k_v))) * (1 / (1 - ksi_pp))
    H_i = kpd_ir * ((h_0 - point_pv.h) + (h_pp - h_1))

    kpd_mech = 0.994
    kpd_eg = 0.98
    G_0 = N_e / (H_i * kpd_mech * kpd_eg * (10 ** 3))
    G_k = (N_e / ((h_k - h_k_v) * kpd_mech * kpd_eg * (10 ** 3))) * ((1 / kpd_ir) - 1)
    return kpd_ir

eta = [Calculate_G0_Gk(p0=p_0, T0=T00, ppp=p_pp, Tpp=T_pp, pk=p_k, Tpv=T_pv, Ne=N_e, deltap0=delta_p_0,
                       deltappp=delta_p_pp, deltap=delta_p) for T00 in T_0]
N_e = 545e6
p_0 = 26.5e6
t_0 = 600
T_0 = t_0+273.15
p_pp = 4.4e6
t_pp = 570
T_pp = t_pp+273.15
p_k = 3.4e3
t_pv = 280
T_pv = t_pv+273.15
delta_p_0 = 0.05*p_0
delta_p_pp = 0.08*p_pp
delta_p = 0.03*p_pp

point_0 = WSP(P=p_0 * 10 ** (-6), T=T_0)
s_0 = point_0.s
h_0 = point_0.h

p_0_ = p_0 - 0.05 * p_0
point_p_0_ = WSP(P=p_0_ * 10 ** (-6), h=h_0)
t_0_ = point_p_0_.T - 273.15
s_0_ = point_p_0_.s

p_1t = p_pp + 0.1 * p_pp
point_1t = WSP(P=p_1t * 10 ** (-6), s=s_0)
t_1t = point_1t.T - 273.15
h_1t = point_1t.h

point_pp = WSP(P=p_pp * 10 ** (-6), T=T_pp)
h_pp = point_pp.h
s_pp = point_pp.s

H_0 = h_0 - h_1t

eta_oi = 0.85
H_i_cvd = H_0 * eta_oi

h_1 = h_0 - H_i_cvd

point_1 = WSP(P=p_1t * 10 ** (-6), h=h_1)
s_1 = point_1.s
T_1 = point_1.T
v_1 = point_1.v

p_pp_ = p_pp - 0.03 * p_pp
point_pp_ = WSP(P=p_pp_ * 10 ** (-6), h=h_pp)
s_pp_ = point_pp_.s
v_pp_ = point_pp_.v

point_kt = WSP(P=p_k * 10 ** (-6), s=s_pp)
T_kt = point_kt.T
h_kt = point_kt.h
v_kt = point_kt.v
s_kt = s_pp

H_0_csdcnd = h_pp - h_kt
eta_oi = 0.85
H_i_csdcnd = H_0_csdcnd * eta_oi
h_k = h_pp - H_i_csdcnd
point_k = WSP(P=p_k * 10 ** (-6), h=h_k)
T_k = point_k.T
s_k = point_k.s
v_k = point_k.v

point_k_v = WSP(P=p_k * 10 ** (-6), x=0)
h_k_v = point_k_v.h
s_k_v = point_k_v.s

eta_oiI = (h_1 - h_0) / (h_1t - h_0)

p_pv = 1.4 * p_0
point_pv = WSP(P=p_pv * 10 ** (-6), T=T_pv)
h_pv = point_pv.h
s_pv = point_pv.s

ksi_pp_oo = 1 - (1 - (T_k * (s_pp - s_k_v)) / ((h_0 - h_1t) + (h_pp - h_k_v))) / (
            1 - (T_k * (s_pp - s_pv)) / ((h_0 - h_1t) + (h_pp - h_pv)))

T_0_ = 374.2 + 273.15

T_ = (point_pv.T - T_k) / (T_0_ - T_k)
if T_ <= 0.636363636:
    ksi1 = -4.9655 * T_ ** 3 + 5.5547 * T_ ** 2 - 1.0808 * T_ + 0.4922
elif 0.636363636 < T_ <= 0.736363636:
    ksi1 = -1.3855 * T_ ** 2 + 2.0774 * T_ + 0.0321
elif 0.736363636 < T_ <= 0.827272727:
    ksi1 = -1.4152 * T_ ** 2 + 2.3287 * T_ - 0.1088
else:
    ksi1 = 0.82

if T_ <= 0.631818182:
    ksi2 = -1.0078 * T_ ** 3 - 0.3296 * T_ ** 2 + 1.7524 * T_ + 0.0714
elif 0.631818182 < T_ <= 0.718181818:
    ksi2 = -2.5821 * T_ ** 2 + 3.689 * T_ - 0.4825
elif 0.718181818 < T_ <= 0.827272727:
    ksi2 = -11.992 * T_ ** 3 + 25.812 * T_ ** 2 - 18.31 * T_ + 5.1453
elif 0.827272727 < T_ <= 0.936363636:
    ksi2 = -11.115 * T_ ** 3 + 27.375 * T_ ** 2 - 22.574 * T_ + 7.1385
else:
    ksi2 = 0.89
ksi = (ksi1 + ksi2) / 2
ksi_pp = ksi * ksi_pp_oo
eta_ir = ((H_i_cvd + H_i_csdcnd) / (H_i_cvd + (h_pp - h_k_v))) * (1 / (1 - ksi_pp))
H_i = eta_ir * ((h_0 - h_pv) + (h_pp - h_1))
eta_m = 0.994
eta_eg = 0.98
G_0 = N_e / (H_i * eta_m * eta_eg * (10 ** 3))
G_k = N_e / ((h_k - h_k_v) * eta_m * eta_eg * (10 ** 3)) * (1 / eta_ir - 1)

itog = pd.DataFrame({
"Температура" : list(range(500, 601, 1)),
"КПД" : (eta)

})
x = (list(range(500, 601, 1)))
y = (eta)

p = figure(
     title='Зависимость КПД от температуры',
     x_axis_label='Температура',
     y_axis_label='КПД')

p.line(x, y, legend_label='Зависимость КПД от температуры', line_width=3)


itog

df = pd.DataFrame({
    "KPD.max" : (max(eta)),
    "Gk" : [G_k],
    "G0" : [G_0]})

df

fig = plt.figure()
point_0 = WSP(P=p_0*1e-6, T=T_0)
p_0_d = p_0 - delta_p_0
point_0_d = WSP(P=p_0_d*1e-6, h=point_0.h)
p_1t = p_pp + delta_p_pp
point_1t = WSP(P=p_1t*1e-6, s=point_0.s)
H_01 = point_0.h - point_1t.h
kpd_oi = 0.85
H_i_cvd = H_01 * kpd_oi
h_1 = point_0.h - H_i_cvd
point_1 = WSP(P=p_1t*1e-6, h=h_1)
point_pp = WSP(P=p_pp*1e-6, T=T_pp)
p_pp_d = p_pp - delta_p_pp
point_pp_d = WSP(P=p_pp_d*1e-6, h=point_pp.h)
point_kt = WSP(P=p_k*1e-6, s=point_pp.s)
H_02 = point_pp.h - point_kt.h
kpd_oi = 0.85
H_i_csd_cnd = H_02 * kpd_oi
h_k = point_pp.h - H_i_csd_cnd
point_k = WSP(P=p_k*1e-6, h=h_k)

s_0 = [point_0.s-0.05,point_0.s,point_0.s+0.05]
h_0 = [WSP(P = p_0*1e-6,s = s_).h for s_ in s_0]
s_1 = [point_0.s-0.05,point_0.s,point_0.s+0.18]
h_1 = [WSP(P=p_1t*1e-6, s = s_).h for s_ in s_1]
s_0_d = [point_0_d.s-0.05, point_0_d.s, point_0_d.s+0.05]
h_0_d = h_0
s_pp = [point_pp.s-0.05,point_pp.s,point_pp.s+0.05]
h_pp = [WSP(P=p_pp*1e-6, s=s_).h for s_ in s_pp]
s_k = [point_pp.s-0.05,point_pp.s,point_pp.s+0.8]
h_k = [WSP(P=p_k*1e-6, s=s_).h for s_ in s_k]
s_pp_d = [point_pp_d.s-0.05,point_pp_d.s,point_pp_d.s+0.05]
h_pp_d = h_pp

plt.plot([point_0.s,point_0.s,point_0_d.s,point_1.s],[point_1t.h,point_0.h,point_0.h,point_1.h],'-or')
plt.plot([point_pp.s,point_pp.s,point_pp_d.s,point_k.s],[point_kt.h,point_pp.h,point_pp.h,point_k.h],'-or')

plt.plot(s_0,h_0)
plt.plot(s_1,h_1)
plt.plot(s_0_d,h_0_d)
plt.plot(s_pp,h_pp)
plt.plot(s_k,h_k)
plt.plot(s_pp_d,h_pp_d)
plt.ylabel('H кДж/кг')
plt.xlabel('S кДж/кг*К')
st.bokeh_chart(p, use_container_width=True)
st.pyplot(fig)