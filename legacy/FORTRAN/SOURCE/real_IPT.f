c############################################################
c IPT solution of the Hubbard model on the real axis at T=0
c by M.J. Rozenberg
c This code is essentially the first code IPT on the real axis
c that was able to go across the Mott-Hubbard transition and
c observe the coexistent solutions. It was written by X.Y. Zhang
c and M.J. Rozenberg, and the results published in X.Y. Zhang,
c M. J. Rozenberg and G. Kotliar, Phys Rev Lett v70, p1666 (1993)
c############################################################

 	parameter(L=2*8192)
	implicit real*8(a-h,o-z)
	double complex one,xi,xs,sq,g,ome
	double complex g0(2*L),tg0(2*L),fg0(-L:L)
	double complex sft2(2*L),sf2(2*L),sft(-L:L),sf(-L:L)
	double precision dns(2*L), dns0(2*L)


c************************************************
c read input data
c f = 	frequency mesh
c d0 = 	half-bandwidth of the semicircular DOS (Bethe lattice)
c u =	Hubbard repulsion

 	read(22,*)f,d0,u
c************************************************

c************************************************
c initialize parameters
 	pi=datan(1.d0)*4.d0
	dtau=2.*pi/f/dfloat(2*L)
	f2=f**2
	one=(1.d0,0.d0)
	xi =(0.d0,1.d0)
	XL=dfloat(L)
	L2=2*L

c set the imet flag for metallic or insulating solution
 	if(u.le.3.d0)imet=1
 	if(u.gt.3.3d0)imet=0
 	if(u.gt.3.d0.and.u.le.3.3d0)then
123	continue
	write(6,*) 'find the metal or the insulator? (1 or 0)'
	read(5,*)imet
	if(imet.ne.0.and.imet.ne.1)then
	write(6,*) 'you have to type 1 for metal or 0 for insulator'
	goto 123
 	endif
 	endif
	
	print*,imet

c nloop= number of sc loops 
	if(imet.eq.0)then
	nloop=30
	else
	nloop=10
	if(u.gt.2.d0) nloop=30
	if(u.gt.3.d0) nloop=100
	if(u.gt.3.2d0) nloop=300
	endif

c nl0= fixes the number of last sc loops to print out
c default is 1
	nl0=nloop-1
	

c output data parameters
	i01=20
	imax1=2000
	i02=30
	imax2=16000

	do i=-L,L
	sf(i)=(0.d0,0.d0)
	enddo

c*******************************************************

c self-consistent loop starts here**********************
	do 100 iloop=1,nloop

	if(iloop.eq.nloop/2)then
	write(6,*) 'I am half way thru...'
	endif


c calculate Go
	do 2 i=L+1,2*L
	om=(float(i)-XL-1.)*f
	sig=1.d0
	ome=om*one

	if(iloop.eq.1)then
c choose a Go seed in the first loop***

	 if(dabs(om).lt.1.d-9)then
	  g0(i)=0.*one
	else
	  g0(i)=1./(ome+float(imet)*d0*xi/2.)
	endif


	else
c calculate a new Go from the last Sigma

	xs=sf(i-L-1)
	sq=cdsqrt((ome-xs)**2-d0**2*one)
	sqim=dimag(sq)
	sqre=dreal(sq)

	if(i.le.L+5)then
c get the first five frequency points on the right branch-cut
		if(sqre.gt.0..and.imet.eq.0)sig=-1.d0
	else
c choose the branch based on continuity criterium 
		benpr=dreal(2./(ome+xs+sq))
	 	benmr=dreal(2./(ome+xs-sq))
	 	benpi=dimag(2./(ome+xs+sq))
		benmi=dimag(2./(ome+xs-sq))
		xp=((benpr+benchr0-2.*benchr)**2
     $    	+(benpi+benchi0-2.*benchi)**2)
		xm=((benmr+benchr0-2.*benchr)**2
     $    	+(benmi+benchi0-2.*benchi)**2)
		if(xp.gt.xm)sig=-1.d0
        	g=(2./(ome+xs+sq))
        	d1=dimag(one/g)-dimag(sf(i-L-1))
        	d2=dreal(one/g)-dreal(sf(i-L-1))
        	d2=d2**2
        	d3=d1**2
        	dn=d1/(d2+d3)
        	if(dn.lt.0)sig=-1.
	endif

c this is the new Go
	g0(i)=2./(ome+xs+sig*sq)

	endif

	benchr=dreal(g0(i))
	benchi=dimag(g0(i))
	benchr0=dreal(g0(i-1))
	benchi0=dimag(g0(i-1))
2	continue


c writes the omega < 0 part of Go for "causal" GF
	do 222 i=1,L-1
	g0(L+1-i)=-g0(L+1+i)
222	continue

c fix the end-point
	g0(1)=g0(2)
c fix the Go(omega=0) equal to 0   (L+1-point) 
	g0(L+1)=0.*one

c*** g0 is Go(omega) 



c-------------------------------------
c FFT from omega to time
        do i=1,2*L
        tg0(i)=g0(i)
        enddo
        call four1(tg0,L2,1)


	ex=-1.
        do 3 i=1,2*L
	ex=-ex
           tg0(i)=ex*f/2./pi*tg0(i)
 3      continue

        do 4 i=1,L
           fg0(i-1)=tg0(i)
 4      continue

        do 5 i=L+1,2*L
           fg0(-2*L-1+i)=tg0(i)
 5      continue

c *** fg0(i) is  Go(t)
c-------------------------------------


c computes the Sigma*************************
c   Sigma(t) =  U^2/4 * Go^3(t)

	do 7 i=-L+1,L-1
	sft(i)=-u**2*fg0(i)**2*fg0(-i)
7	continue
	sft(-L)=-u**2*fg0(-L)**2*fg0(L-1)

        do 24 i=-L,L-1
	   sft2(i+L+1)=sft(i)
 24      continue


c********************************************

c-------------------------------------
c FFT from time to omega 

        do i=1,2*L
        sf2(i)=sft2(i)
        enddo
        call four1(sf2,L2,-1)


	ex=-1.
	do 8 i=1,2*L
	ex=-ex
	sf2(i)=-ex*dtau*sf2(i)
8	continue

        do 34 i=1,L
           sf(i-1)=sf2(i)
34      continue

        do 35 i=L+1,2*L
           sf(-2*L-1+i)=sf2(i)
35      continue


c*** sf is Sigma(omega)
c-------------------------------------


c print output data********************************************

        if(iloop.eq.1.and.imet.eq.1)then
c get the first loop Im[G(w)]  (impurity problem, no sc)*********
        do 110 i=1,2*L
        if(i.eq.L+1)then
c computes the value of Im[G(0)]
        dns0(i)=2.*float(imet)/d0
        else
c computes the value of Im[G(w)]
        d1=dimag(one/g0(i))-dimag(sf(i-L-1))
        d2=dreal(one/g0(i))-dreal(sf(i-L-1))
        d2=d2**2
        d3=d1**2
        dns0(i)=d1/(d2+d3)
        endif
110     continue
	endif


	if(iloop.gt.nl0)then
c get the self-consistent Im[G(w)]****************
	do 111 i=1,2*L
	if(i.eq.L+1)then
c computes the value of Im[G(0)]
	dns(i)=2.d0*dfloat(imet)/d0
	else
c computes the value of Im[G(w)]
	d1=dimag(one/g0(i))-dimag(sf(i-L-1))
	d2=dreal(one/g0(i))-dreal(sf(i-L-1))
	d2=d2**2
	d3=d1**2
	dns(i)=d1/(d2+d3)
	endif
111	continue



c print output files
c fort.25 is Im[G(w)]   (the DOS)
c fort.23 is Re[Sigma(w)]
c fort.24 is Im[Sigma(w)]
c fort.26 is Im[G(w)]  (impurity problem, no self-consistency, imet=1)

	do  i=imax2,imax1+1,-i02
        write(26,*)-real(f*float(i)),(dns0(i+L+1))
        write(25,*)-real(f*float(i)),(dns(i+L+1))
 	write(24,*)-real(f*float(i)),dabs(dimag(sf(-i)))
        write(23,*)-real(f*float(i)),(dreal(sf(-i)))
c	write(27,*)-real(f*float(i)),dreal(dreal(g0(i+L+1)))
c	write(28,*)-real(f*float(i)),dreal(dimag(g0(i+L+1)))
	enddo

        do i=imax1,0,-i01
        write(26,*)-real(f*float(i)),(dns0(i+L+1))
        write(25,*)-real(f*float(i)),(dns(i+L+1))
 	write(24,*)-real(f*float(i)),dabs(dimag(sf(-i)))
        write(23,*)-real(f*float(i)),(dreal(sf(-i)))
c	write(27,*)-real(f*float(i)),dreal(dreal(g0(i+L+1)))
c	write(28,*)-real(f*float(i)),dreal(dimag(g0(i+L+1)))
	enddo

	do 118 i=1,imax1,i01
        write(26,*)real(f*float(i)),(dns0(i+L+1))
	write(25,*)real(f*float(i)),(dns(i+L+1))
 	write(24,*)real(f*float(i)),dabs(dimag(sf(i)))
 	write(23,*)real(f*float(i)),(dreal(sf(i)))
c	write(27,*)real(f*float(i)),dreal(dreal(g0(i+L+1)))
c	write(28,*)real(f*float(i)),dreal(dimag(g0(i+L+1)))
118	continue

	do 116 i=imax1+1,imax2,i02
        write(26,*)real(f*float(i)),(dns0(i+L+1))
	write(25,*)real(f*float(i)),(dns(i+L+1))
 	write(24,*)real(f*float(i)),dabs(dimag(sf(i)))
 	write(23,*)real(f*float(i)),(dreal(sf(i)))
c	write(27,*)real(f*float(i)),dreal(dreal(g0(i+L+1)))
c	write(28,*)real(f*float(i)),dreal(dimag(g0(i+L+1)))
116	continue

 	write(23,*)'    '
 	write(24,*)'    '
	write(25,*)'    '
 	write(26,*)'    '
c	write(27,*)'    '
c	write(28,*)'    '
	endif
c*************************************************************


c go back to the self consistency loop
100	continue	


	stop
	end

c########################################################
      SUBROUTINE four1(data,nn,isign)
        implicit real*8(a-h,o-z)
        double precision data(2*nn)
c     INTEGER isign,nn
c     REAL data(2*nn)
c     INTEGER i,istep,j,m,mmax,n
c     REAL tempi,tempr
c     DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END

c########################################################
c########################################################


