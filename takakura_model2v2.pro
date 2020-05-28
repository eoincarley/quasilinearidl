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

function sfun, vindex, tindex, xindex, $
		Ftot=Ftot, Wtot=Wtot, Fthermal=Fthermal, Fs=Fs, Wthermal=Wthermal, $
		const0=const0, const1=const1

	dFtot = Ftot[vindex+1, tindex, xindex] - Ftot[vindex, tindex, xindex]
	dFtherm = Fthermal[vindex+1]-Fthermal[vindex]
	valS = const0*( dFtot*Wtot[vindex, tindex, xindex] - dFtherm*Wthermal[vindex] ) $
		+ const1*Fs[vindex, tindex, xindex]
	
	;if ~finite(valS) then stop
	return, valS
end

pro takakura_model2v2, postscript=postscript


root = '/Users/eoincarley/python/quasilinear'


;		Constants

kb = 1.38d-16	; egs/K
e = 4.803d-10	; StatCoulomb
me = 9.109d-28	; g
omega = 1.0d8  	; Hz
ndensity = 1.235d8    ; cm^-3 (np.pi*me*omega**2)/e**2
Te = 1d6 		; K
psi = 6d0 		; Electron spectral index
vd = 2d9 		; cm/s
Ns = 4e9		; cm^-2
dreal = 1e9 	; cm
L = 1.6002e11	; cm
vthermal = sqrt(2d0*kb*Te/me)  ; cm/s
vc = vthermal/sqrt(2.0) 		; cm/s
vc = vd 		; cm/s
c = 2.9e10 		; cm/s
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
F = dblarr( n_elements(V), n_elements(T), n_elements(X) )
Fp = dblarr( n_elements(V), n_elements(T), n_elements(X) )
Ff = dblarr( n_elements(V), n_elements(T), n_elements(X) )

for m=0, n_elements(V)-1 do begin
	for l=0, n_elements(X)-1 do begin
		for n=0, n_elements(T)-1 do begin
			Ff[m, n, l] = exp( -( (X[l]-V[m]*T[n])/D )^2.0 )*G*V[m]^(-psi)
		endfor
	endfor
endfor	

Fs[1:nvsteps-1, *, *] = Ff[1:nvsteps-1, *, *] + Fp[1:nvsteps-1, *, *]
for i=0, n_elements(X)-1 do begin
	for j=0, n_elements(T)-1 do begin
		F[*, j, i] = Fs[*, j, i] + FT[*]
	endfor
endfor	

S = dblarr( n_elements(V), n_elements(T), n_elements(X) )
W = dblarr( n_elements(V), n_elements(T), n_elements(X) )

wthermal = (0.25*me*omega^2)*(vthermal/vc)^2 *(1d0-(vc/vreal)^2)
wthermal[where(vreal lt vc)] = 0.0
WT = wthermal/(me*omega^2)
for i=0, n_elements(X)-1 do W[*, 0, i] = WT[*]


Gamma1 = dblarr( n_elements(V) )
Gamma2 = dblarr( n_elements(V) )
Phi = dblarr( n_elements(V) )
decay = dblarr( n_elements(V) )

c2a = 2d0*delt*V
c2b = delt*V
c3 = delv*V
c3b = 2d0/(delv*V)
c4 = delt*V*alog(V)
c6 = WT*(delt/(delv*Omega_p)) * V^2.0
c6a = delt*Cf*WT
c7 = (delt/delv) * V
tindex0 = round(11)
tindex1 = round(10)

sconst0 = (1.0/(delv*V))
sconst1 = (Omega_p^2/V^2)*(alog(V) + 2.0*alog(Lambda))

loadct, 1
window, 0, xs=800, ys=1000
window, 1, xs=400, ys=400
mi1 = where(V lt 1.2 and V gt 1.0)
nstart=1
for l=1, n_elements(X)-1 do begin

	for n=nstart, n_elements(T)-2 do begin

		; Change X-grid spacing after a certain number of iterations.
		if l le 100 then delx = delx0
		if l gt 100 and l le 150 then delx = delx1	
		if l gt 150 then delx = delx2	

		;------------------------------------------;
		;
		;		Electron beam distribution
		;
		for m=1, n_elements(V)-2 do begin
		
			S1 = sfun(m, n, l-1, Ftot=F, Wtot=W, Fthermal=FT, Fs=Fs, Wthermal=WT, const0 = sconst0[m], const1 = sconst1[m])
			;if ~finite(S1) then stop
			S0 = sfun(m-1, n, l-1, Ftot=F, Wtot=W, Fthermal=FT, Fs=Fs, Wthermal=WT, const0 = sconst0[m], const1 = sconst1[m])
			;if ~finite(S1) then stop
			delS = S1-S0

			if l eq 1 or X[l]/D>6.1 then begin
				Fp[m, n, l] = Fp[m, n, l-1] + $
						 ( delx/c2a[m] ) * (  Fp[m, n-1, l-1] - Fp[m, n+1, l-1] ) + $
						 ( delx/c3[m] ) * delS
			endif else begin			 
				Fp[m, n, l] = Fp[m, n, l-2] + $
							 ( delx/c2b[m] ) * (  Fp[m, n-1, l-1] - Fp[m, n+1, l-1] ) + $
							 ( c3b[m]*delx ) * delS
			endelse				
		
			Fs[m, n, l] = Ff[m, n, l] + Fp[m, n, l]
			F[m, n, l] = Fs[m, n, l] + FT[m]
			
		endfor

		;if X[l]/D ge 3.6 and n eq 154 then stop
		
		;------------------------------------------;
		;
		;	   Langmuir wave energy distribution
		;
		for m=1, n_elements(V)-2 do begin

			Phi[m] = c4[m]*Fs[m, n, l] - c6[m]*(FT[m+1]-FT[m]) + c6a[m]
			Gamma1[m] = c7[m]*( F[m+1, n, l] - F[m, n, l]  ) - delt*Cf
			Gamma2[m] = Phi[m]/W[m, n-1, l]
		
			if abs(Gamma1[m]) gt (Gamma2[m] + 2e-1) then $
				W[m, n, l] = W[m, n-1, l]*exp(Gamma1[m] + Gamma2[m])

			if abs(Gamma1[m]) le (Gamma2[m] + 2e-1) then $
				W[m, n, l] = (1d0 + Gamma1[m])*W[m, n-1, l] + Phi[m]

			; Check for the decay phase at each velocity
			if m gt 2 and n gt tindex0 then begin
				dW = deriv(T[n-tindex1:n], W[m, n-tindex1:n, l])
				negdW = where(dW lt 0.0)
				if n_elements(negdW) eq tindex0 or decay[m] eq 1 then begin
					W[m, n, l] = -Phi[m]/Gamma1[m]
					;print, 'Decay'
					decay[m]=1
				endif else begin
					;print, 'Rise'
				endelse	
			endif 
	
		endfor	

		W[mi1, n, l] = -Phi[mi1]/Gamma1[mi1]


		nplot = 200/tres
		if X[l]/D ge 4.0 and n lt nplot then begin
			wset, 0
			plot_image, congrid( sigrange(reform(Fs[*, 0:nplot, l-1])), ntsteps, ntsteps), pos = [0.1, 0.7, 0.5, 0.95], /noerase, title='Fs'
			plot_image, congrid( bytscl(reform(Fp[*, 0:nplot, l]), -50, 50), ntsteps, ntsteps), pos = [0.1, 0.4, 0.5, 0.65], /noerase, title='Fp'
			plot_image, congrid( sigrange(reform(alog10(F[*, 0:nplot, l]))), ntsteps, ntsteps), pos = [0.1, 0.1, 0.5, 0.35], /noerase, title='F', xtitle='Velocity', ytitle='time'	
			plot_image, congrid( sigrange(reform(alog10(W[*, 0:nplot, l]) )), ntsteps, ntsteps), pos = [0.55, 0.7, 0.95, 0.95], /noerase, title='W'

			set_line_color
			wset, 1
			plot, V, F[*, n, l], pos = [0.15, 0.55, 0.95, 0.95], yr=[1e-6, 1e4], /ylog, title=X[l]/D, color=1
			oplot, V, Fp[*, n, l], color=5
			plot, V, W[*, n, l], pos = [0.15, 0.1, 0.95, 0.5], yr=[-10, 10], title=X[l]/D, /noerase

			loadct, 1, /silent
			if n eq 40 then stop
			
		endif
		;if X[l]/D ge 3.6 then nstart=20
		;if n eq 1 then stop
	endfor
	decay[*]=0.0
	progress_percent, l, 1, n_elements(X)-1	
	
	
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
	in = (where(X/D ge 7.0))[0]

	vin = where(V eq [1.4, 1.6, 1.8, 2.0]) 
	cols = interpol([0,250], n_elements(vin))
	plot, T, Fs[vin[0], *, xin], /ylog, yr=[10, 1e5], /xs, $
		pos=[0.2, 0.5, 0.9, 0.9], /normal, /noerase, /nodata, thick=3

	plot, T, Fs[vin[0], *, xin], /ylog, yr=[10, 1e5], /xs, $
		pos=[0.2, 0.5, 0.9, 0.9], /normal, /noerase, thick=3, $
		xtitle = 'Time', ytitle='Fs'
	loadct, 33
	for i=0, n_elements(vin)-1 do begin

		loadct, 33, /silent
		oplot, T, Fs[vin[0], *, xin], color=cols[i], thick=3

	endfor

	

if keyword_set(postscript) then device, /close
	set_plot, 'x'

stop
END