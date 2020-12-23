function xout = czt_2d( xin , scalx , scaly , idir )
%
%  Aufruf:
%  xout = czt_two(xin,scalx,scaly,idir);
%
%  Beschreibung:
%  Berechnung der zweidimensionalen Chirp-z-Transformation mittels Faltung
%  mit Skalierungsfaktoren scalx bzw. scaly f黵 ein zweidimensionales Feld
%  xin. Beide Koordinatenschnitte werden unabh鋘gig behandelt. 
%  Das Feld kann rechteckig sein mit Dimensionen npx / npy.
%  Basis ist die Routine czt aus der Signal-Toolbox von Matlab.
%  F黵 scalx/scaly = 1 entsteht die FFT. Die Frequenzachsen werden um die Faktoren
%  scalx/scaly gespreizt.
%  Es mu� scalx/scaly => 1 gelten.
%  Die Index-Shiftfunktion ist bereits enthalten, die Frequenz s = 0 liegt beim Index np/2+1.
%  Es gilt f黵 die Frequenzrasterweite :    ds = 1 / ( npx * dx * scalx )  (y analog)
%

%% English 
% Calculation of two-dimensional chirp z-transform means convolution
% With scaling factors scalx or scaly for a two-dimensional field
% Xin. Both coordinates cuts are treated independently.
% The field can be rectangular with dimensions npx / NPY.
% Base is the routine czt from the signal toolbox of Matlab.
% For scalx / scaly = 1 arises the FFT. The frequency axes are the factors
% Scalx / spread scaly.
% It must scalx / scaly apply => 1.
% The index shift function is included, the frequency of s = 0 is the index np / 2 +1.
% It applies to the frequency screen width: ds = 1 / (dx * scalx npx *) (y analog)
%%

%  Version:
%  23.06.06   Herbert Gross  7.0  躡erarbeitete Version 
%
%  Input:   xin(npx,npy)  : Inputfeld, zweidimensional
%           scalx         : Lupenfaktor in x-Richtung
%           scaly         : Lupenfaktor in y-Richtung
%           idir          : Steuerparameter
%                           idir = 0 : Transformation in x und y
%                           idir = 1 : Transformation nur in x-Richtung
%                           idir = 2 : Transformation nur in y-Richtung
%
%  Output:  xout(npx,npy) : Outputfeld, zweidimensional
%
%  Abh鋘gigkeiten : keine
%
%  Referenzen:
%  Rabiner, Schafer, Rader, The Chirp z-Transform Algorithm, IEEE Trans AU 17(1969) p. 86
%
%  Beispiel: Transformation eines Gaussprofils, verschiedene Punktzahlen und 
%  Spreizungsfaktoren
%
%  npx = 256  ; npy = 128 ;
%  scalx = 8 ;  scaly = 4 ; idir = 0 ;
%  a = 1; xmax = 2.;xmin = -xmax;
%  dx = ( xmax-xmin ) / (npx-1); dy = ( xmax-xmin ) / (npy-1);
%  for j=1:npx ;   x(j) = xmin + (j-1)*dx; end
%  for k=1:npy ;   y(k) = xmin + (k-1)*dy ; end  
%  [yp xp] = meshgrid(y,x); rp = sqrt(xp.^2+yp.^2);indin = find( rp < a );
%  xin = zeros(npx,npy,1); xin(indin) = 1 ;
%  xout = czt_2d( xin , scalx , scaly , idir );
%  figure(1);
%  pcolor(abs(xout));shading flat
%___________________________________________________________________________________
%
[ npx , npy ] = size( xin );
xout = zeros(npx,npy,1);
%
%  x-Richtung, 1. Index
%
if idir == 0 || idir == 1
%
   fak = zeros(1,npx);
   n2 = npx/2+1;
   nn = (0:(npx-1))';
   kk = ( (-npx+1):(npx-1) ).';
   kk2 = (kk .^ 2) ./ 2;
%
   w = exp( -1i*2*pi/(npx*scalx));
   a = exp( -1i*pi/scalx );
   ww = w .^ (kk2);   
   aa = a .^ ( -nn );
   aa = aa .* ww(npx+nn);
   fv = fft( 1 ./ ww(1:(2*npx-1)), 2*npx );   
   for j=1:npx ; fak(j) =  exp(pi*1i*(j-n2)/scalx );end
   fak = ww( npx:(2*npx-1) ) .* fak.' ;
%
   for k=1:npy
      y = xin(:,k) .* aa;
      g = ifft( fft(  y , 2*npx ) .* fv );
      xout(:,k) = ( g( npx:(2*npx-1) ) .* fak ) ;
   end
%
else
    xout = xin ;
%
end
%
%  y-Richtung, 2. Index
%
if idir == 0 || idir == 2
%
   clear fak ;
   fak = zeros(1,npy);
   n2 = npy/2+1;
   nn = (0:(npy-1))';
   kk = ( (-npy+1):(npy-1) ).';
   kk2 = (kk .^ 2) ./ 2;
   w = exp( -1i*2*pi/(npy*scaly));
   a = exp( -1i*pi/scaly );
   ww = w .^ (kk2);   
   aa = a .^ ( -nn );
   aa = aa .* ww(npy+nn);
   fv = fft( 1 ./ ww(1:(2*npy-1)), 2*npy );   
   for j=1:npy ; fak(j) =  exp(pi*1i*(j-n2)/scaly );end
   fak = ww( npy:(2*npy-1) ) .* fak.' ;
%
   for j=1:npx
      y = xout(j,:).' .* aa;
      g = ifft( fft(  y , 2*npy ) .* fv );
      xout(j,:) = ( g( npy:(2*npy-1) ) .* fak ).' ;
   end
%
end
%
xout = conj(xout);