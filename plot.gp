reset

set key top left reverse Left
set grid
set tics out nomirr

set term pdf enh
set out "calib.pdf"



# 0.5 ns	ns	0.474
# 1 ns	    ns	1.179
# 2 ns	    ns	1.995
# 5 ns	    ns	5.918
# 10 ns	    ns	10.947
# 20 ns	    ns	20.694
# 50 ns	    ns	49.100
# 100 ns	ns	98.526
# 200 ns	ns	191.025
# 500 ns	ns	507.852
# 1 us	    us	1003.027
# 2 us	    us	1990.008
# 5 us	    us	5076.382
# 10 us	    us	9919.118
# 20 us	    us	19816.223
# 50 us	    us	48585.712
# 100 us	us	97750.030
# 200 us	us	197578.876
# 500 us	us	499607.410
# 1 ms	    ms	1001235.518



gauss(x)=gauss_a+gauss_b*exp(-(0.5*((x-gauss_center)/gauss_sigma)**2))
fwhm(sig)=2*sqrt(2*log(2))*sig
gauss_a = 1774.95379976011
gauss_b = 9449.75394237946
gauss_sigma = 3.34134367520707
gauss_center = 415.705984497551

fileslit='slit-calib.txt'
fit gauss(x) fileslit i 1 via gauss_sigma,gauss_center,gauss_b,gauss_a

FWHM=fwhm(gauss_sigma)

# we will use this as error 
ratio=1024./FWHM/2

set title sprintf("FWHM = %f",FWHM)

set xl "pixel"
set yl "counts"

plot [gauss_center-50:gauss_center+50] \
fileslit i 1 w l title sprintf("%s",fileslit), \
gauss(x) title "gaussian fit"

print "FWHM=",FWHM



MyFileSlit="slit-before61.txt"
fit gauss(x) MyFileSlit i 1 via gauss_sigma,gauss_center,gauss_b,gauss_a
my_FWHM=fwhm(gauss_sigma)

set title sprintf("FWHM = %f",fwhm(gauss_sigma))
plot [gauss_center-50:gauss_center+50] \
MyFileSlit i 1 w l title sprintf("%s",MyFileSlit), \
gauss(x) title "gaussian fit"

print "my FWHM=",my_FWHM


unset title


linfit(x)= lamp_a*x+lamp_b

fit linfit(x) 'results' u 5:3:($4*ratio) yerror via lamp_a,lamp_b 

set title sprintf("counts=%fx + %f",lamp_a,lamp_b)

set xl "time window [ns]"
set yl "graybody counts (raw)"
set xr [0:1.1e6]
plot linfit(x) lw 2 notitle , \
'results' u 5:3:($4*ratio) w yerr pt 0 notitle


set log xy

set xr [1:1.1e6]
set yl "graybody counts (without offset)"

plot (linfit(x)-lamp_b) lw 2 notitle ,\
'results' u 5:($3-lamp_b):($4*ratio) w yerr pt 1 title sprintf("x%f(+/-%f)",lamp_a,lamp_a_err)

unset log xy


filter=450.0
filter_err=2

unset title




alfa=1.439E-2/(filter*1E-9) # hc /(kb lambda)

Tlamp=2764 #

emissivity=0.457 #rapport emissivity(450nm)/emissivity(650nm)

TBB=alfa/(log(1.0+(exp(alfa/Tlamp)-1.0)/emissivity))


time_window=20.694

k=(exp(alfa/TBB)-1)*lamp_a/FWHM


T0=alfa/11604.5

T(counts,slit,time,rifl)=T0/(log(1.0+(1.0-rifl)*slit*time_window*k/counts))

set xr [0:20000]

titolo=sprintf("T(slit,counts,expTime,R)=T0/(log(1+(1.0-R)*A/counts))\n T0=%.4feV=%.0fK; A=k*slit*expTime\n k=%.4f; slit=%.2fpx; time=%.2fns",T0,T0*11604.5,k,my_FWHM,time_window)


set title titolo.sprintf(" A=%.0f",k*my_FWHM*time_window)

set xl "Counts"
set yl "Temperature eV"


plot \
for [Refl = 0:50:10] T(x,my_FWHM,time_window,Refl/100.0)  lw 3 title sprintf("R=%d\%",Refl)

set xl "Counts (R=0)"
plot \
for [slit = 5:10:1] T(x,slit,time_window,0)  lw 3 title sprintf("slit=%dpx A=%.0f",slit,k*slit*time_window)

set xr[0:1000]
plot \
T(x,my_FWHM,time_window,0.0)*11604.5  lw 3 title sprintf("R=%d\%",0.0)

set xr[0:1000]
plot \
for [new_FWHM = 5:20:5] T(x,new_FWHM,time_window,0.0)*11604.5  lw 3 title sprintf("slit FWHM=%dpx",new_FWHM)




show var

