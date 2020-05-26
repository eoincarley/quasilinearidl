pro setup_ps, name, xsize, ysize

    set_plot,'ps'
    !p.font=0
    !p.charsize=1.0
    device, filename = name, $
        ;/decomposed, $
        /color, $
        /helvetica, $
        /inches, $
        xsize=xsize, $;xsize/100, $
        ysize=ysize, $;xsize/100, $
        /encapsulate, $
        bits_per_pixel=64;, $
       ; yoffset=5

end

pro takakura_model2, postscript=postscript


root = '/Users/eoincarley/python/quasilinear'


;		Constants

kb = 1.38d-16	; egs/K
e = 4.8d-10		; StatCoulomb
me = 9.109d-28	; g
omega = 1.0d8  	; Hz
ndensity = 1.235d8    ; cm^-3 (np.pi*me*omega**2)/e**2
Te = 1d6 		; K
psi = 6d0
vd = 2d9 		; cm/s
Ns = 4e10		; cm^-2
dreal = 1d9 		; cm
L = 1.6002d11	; cm
vthermal = sqrt(2d0*kb*Te/me)  ; cm/s
vc = vthermal/sqrt(2) 		; cm/s
c = 2.9d10 		; cm/s
VT = vthermal/vd
Omega_p = 1.0
Lambda = 3.0*sqrt( (kb^3*Te^3)/(!pi*ndensity) )/(2.0*e^3)
Cf = 0.777*alog(Lambda)*Omega_p^2/VT^3
Sc = 0.0
;--------------------------#
;	Normalised variables
;	and initial conditions
;
tres = 1.0 ;1d-0
xres = 1.0 ;1d-0
vres = 1.0 ;*tres

print, 'V and T res ratio: ' + string(vres/tres)

D = dreal/L
delv = 1d-1*vres
delt = (3d-2)*D*tres
delx0 = (2d-2)*D*xres
delx1 = (1d-2)*D*xres
delx2 = (5d-3)*D*xres

nvsteps = 15.0/vres
ntsteps = 370.0/tres
nxsteps0 = 100.0/xres
nxsteps1 = 50.0/xres
nxsteps2 = 300.0/xres
;------------------------------;
;
;			X, V, T
;
X0 = dindgen(nxsteps0)*delx0 + 3.4*D
X1 = dindgen(nxsteps1)*delx1 + 5.4*D
X2 = dindgen(nxsteps2)*delx2 + 5.9*D
X = [X0, X1, X2]
V = dindgen(nvsteps)*delv + 1.0
vreal = V*vd
T = dindgen(ntsteps)*delt


;------------------------------;
;
;	 All the F distributions
;
G = ((psi-1.0)*Ns)/(dreal*2.0*!pi^(3./2.))*(vd/omega)^3
Gth = ndensity/(VT*2d0*!pi^(3./2.))*(vd/omega)^3
FT = Gth*exp(-1d0*(V/VT)^2.0)
Fs = dblarr( n_elements(V), n_elements(T), n_elements(X) )

for m=1, n_elements(V)-1 do begin
	for l=0, n_elements(X)-1 do begin
		Fs[m, 1, l] = exp(-(X[l]/D)^2)*G*V[m]^(-psi)
	endfor
endfor	

Ff = dblarr( n_elements(V), n_elements(T), n_elements(X) )
for m=0, n_elements(V)-1 do begin
	for l=0, n_elements(X)-1 do begin
		for n=0, n_elements(T)-1 do begin
			Ff[m, n, l] = exp( -( (X[l]-V[m]*T[n])/D )^2.0 )*G*V[m]^(-psi)
		endfor
	endfor
endfor	

F = dblarr( n_elements(V), n_elements(T), n_elements(X) )
Fp = dblarr( n_elements(V), n_elements(T), n_elements(X) )

Fs = Ff + Fp
for i=0, n_elements(X)-1 do begin
	for j=0, n_elements(T)-1 do begin
		F[*, j, i] = Fs[*, j, i] + FT[*]
	endfor
endfor	




S = dblarr( n_elements(V), n_elements(T), n_elements(X) )
W = dblarr( n_elements(V), n_elements(T), n_elements(X) )

wthermal = (0.25*me*omega^2)*(vthermal/vc)^2 *(1d0-(vc/vreal)^2)
;wthermal[where(vreal lt vc)] = 0.0
WT = wthermal/(me*omega^2)
for i=0, n_elements(X)-1 do W[*, 0, i] = WT[*]

Gamma1 = dblarr( n_elements(V) )
Gamma2 = dblarr( n_elements(V) )
Phi = dblarr( n_elements(V) )

dFT = FT[2:-1] - FT[1:-2]
dFTWT = dFT*WT[1:-2]
Vind = V[1:-2]
c0 = ( 1d0/(delv*Vind) )
c1 = (Omega_p^2/Vind^2)*(alog(Vind) + 0.0)
c2 = delt*V[2:-2]
c3 = delv*V[2:-2]
c4 = delt*Vind*Omega_p*alog(Lambda)
c5 = delt*((WT[1:-2]*Vind^2)/(delv*Omega_p) )* dFT
c6 = (delt/delv) * (Vind^2.0/Omega_p^2.0)

loadct, 1
window, 0, xs=500, ys=700
stop
for l=1, n_elements(X)-1 do begin
	maxG=0
	for n=1, n_elements(T)-2 do begin

		if l le 100 then delx = delx0
		if l gt 100 and l le 150 then delx = delx1	
		if l gt 150 then delx = delx2	


		;--------------------------------------;
		;
		;			Beam evolution
		;
		dF = F[2:-1, n, l-1] - F[1:-2, n, l-1]
		S[1:-2, n, l-1] = c0*( dF*W[1:-2, n, l-1] - dFTWT ) + c1*Fs[1:-2, n, l-1]
		delS = S[2:-2, n, l-1] - S[1:-3, n, l-1]

		if l eq 1 or X[l]/D>6.1 then begin
			Fp[2:-2, n, l] = Fp[2:-2, n, l-1] + $
					 ( delx/(2d0*c2) ) * (  Fp[2:-2, n-1, l-1] - Fp[2:-2, n+1, l-1] ) + $
					 ( delx/c3 ) * delS
		endif else begin			 
			Fp[2:-2, n, l] = Fp[2:-2, n, l-2] + $
						 ( delx/c2 ) * (  Fp[2:-2, n-1, l-1] - Fp[2:-2, n+1, l-1] ) + $
						 ( 2d0*delx/c3 ) * delS
		endelse				

		Fs[1:-2, n, l] = Ff[1:-2, n, l] + Fp[1:-2, n, l]
		F[1:-2, n, l] = Fs[1:-2, n, l] + FT[1:-2]

		;--------------------------------------;
		;
		;	Calculte Langmuir wave energy
		;
		Phi[1:-2] = c4*Fs[1:-2, n, l] - c5 + delt*Cf*WT[1:-2]
		Gamma1[1:-2] = c6*( F[2:-1, n, l] - F[1:-2, n, l] ) - delt*Cf
		Gamma2[1:-2] = Phi[1:-2]/W[1:-2, n-1, l]

		mi1 = where( abs(Gamma1) gt (Gamma2 + 2e-1) )
		if mi1[0] gt -1 then $
			W[mi1, n, l] = W[mi1, n-1, l]*exp(Gamma1[mi1] + Gamma2[mi1]) 

		mi2 = where( abs(Gamma1) le (Gamma2 + 2e-1) )
		if mi2[0] gt -1 then $	
			W[mi2, n, l] = (1d0 + Gamma1[mi2])*W[mi2, n-1, l] + Phi[mi2]


		mi1 = where(V lt 1.2)
		if n_elements(mi1[0]) gt 0 then $
			W[mi1, n, l] = -Phi[mi1]/Gamma1[mi1]

		; Check for the decay phase at each velocity
		;for m=3, n_elements(V)-1 do begin
		;	decay_phase = 0
		;	if n gt round(11/tres) then begin
		;		dW = deriv(T[n-round(10/tres):n], W[m, n-round(10/tres):n, l])
		;		negdW = where(dW lt -1.0)
		;		if n_elements(negdW) eq round(11/tres) then decay_phase=1
		;	endif	
;
;			if decay_phase eq 1 and W[m, n, l] lt 0.1 then begin
;				if W[m, n, l] lt 0.1 then $
;					W[m, n, l] = -Phi[m]/Gamma1[m]
;			endif 
;		endfor	


		if X[l]/D ge 0.0 and n lt 200 then begin
			wset, 0
			;plot_image, congrid( sigrange(reform(alog10(Ff[*, *, l-1]))), 100, 100), pos = [0.1, 0.75, 0.5, 0.95], /noerase
			plot_image, congrid( bytscl(reform(Fp[*, *, l]), -10, 10), ntsteps, ntsteps), pos = [0.1, 0.5, 0.5, 0.75], /noerase
			;plot_image, congrid( sigrange(reform(Fs[*, *, l-1])), 100, 100), pos = [0.1, 0.35, 0.5, 0.55], /noerase
			plot_image, congrid( sigrange(reform(alog10(F[*, *, l]))), ntsteps, ntsteps), pos = [0.1, 0.1, 0.5, 0.35], /noerase, xtitle='Velocity', ytitle='time'	
			plot_image, congrid( sigrange(reform(alog10(W[*, *, l]) )), ntsteps, ntsteps), pos = [0.55, 0.75, 0.95, 0.95], /noerase
			
			set_line_color
			wset, 1
			;plot, Vind, F[1:-2, n, l], yr=[1, 1e3], title = string(maxG1)+' '+string(X[l]/D), pos = [0.15, 0.55, 0.95, 0.95]
			;maxG1 = max([maxG1, max(Gamma1, /nan)])
			;plot, Vind, Gamma1, yr=[-0.6, 0.6], pos = [0.15, 0.15, 0.95, 0.5], /noerase
			;oplot, [1.0, 2.4], [maxG1, maxG1], color=3
			;oplot, [1.0, 2.4], [0.0, 0.0], linestyle=1

			plot, Vind, Fp[*, n, l], pos = [0.15, 0.15, 0.95, 0.95], yr=[-2,2]

			loadct, 1, /silent
			;wait, 0.5
			;stop
			
		endif	

	endfor
	
	progress_percent,l,1,n_elements(X)-1	
	
endfor	

loadct, 0

if keyword_set(postscript) then begin
	setup_ps, '~/Desktop/quasilinear_diffusion.eps', 6, 9
endif else begin
	window, 0, xs=500, ys=700
endelse
	;----------------------------------;
	;
	;	  Particle distributions
	;
	xin = (where(X/D ge 7.0))[0]
	;tin = [70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 140, 170]
	;tin = [95, 105, 115, 125, 135, 145]/tres
	tin = [115, 125, 133, 142, 149, 157, 170]/tres
	cols = interpol([0,250], n_elements(tin))
	plot, V, Ff[*, tin[0], xin], /ylog, yr=[10, 1e4], xr=[1.1, 2.0], /xs, $
		pos=[0.2, 0.5, 0.9, 0.9], /normal, /noerase, /nodata, thick=3

	plot, V, Ff[*, tin[0], xin], /ylog, yr=[10, 1e4], xr=[1.1, 2.0], /xs, $
		pos=[0.2, 0.5, 0.9, 0.9], /normal, /noerase, $
		xtitle = ' ', ytitle='Number density (N)'
	loadct, 33
	for i=0, n_elements(tin)-1 do begin
		
		loadct, 0, /silent
		oplot, V, Ff[*, tin[i], xin], linestyle=2, color=100, thick=1

		loadct, 33, /silent
		oplot, V, F[*, tin[i], xin], color=cols[i], thick=3
		oplot, V, Fp[*, tin[i], xin], color=cols[i], linestyle=1, thick=3
		
	endfor
	loadct, 0
	oplot, V, FT, linestyle=1, thick=3, color=0


	;----------------------------------;
	;
	;	  Langmuir distributions
	;
	;tin = [85, 95, 100, 105, 110, 115, 120, 130, 135, 139, 145, 170]/tres
	tin = [115, 125, 133, 142, 149, 157, 161, 170]/tres
	cols = interpol([0,250], n_elements(tin))

	plot, V, W[*, tin[0], xin], /ylog, yr=[1e-1, 1e2], xr=[1.1, 2.0], $
		/xs, pos=[0.2, 0.1, 0.9, 0.5], /normal, /noerase, $
		xtitle = 'Velocity (V)', ytitle='Langmuir energy (W)'
	loadct, 33
	for i=0, n_elements(tin)-1 do begin
		
		oplot, V, W[*, tin[i], xin], color=cols[i], thick=3
		
	endfor
	loadct, 0
	oplot, V, WT, linestyle=1, thick=3

if keyword_set(postscript) then device, /close
	set_plot, 'x'

stop
END