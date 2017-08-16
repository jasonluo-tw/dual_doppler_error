       Program plot_dual_doppler_outdata
       integer ixmax,iymax,domain
       character filein*100       
C------------------------------------------------------------

       write(*,*)'input CEDREAD outdata ='
c       read(*,'(A100)')filein
       filein='case0257'
       open(10,file=filein,status='old')
       read(10,*)ixmax,iymax,level,ifield
       write(*,*)'Choose domain'
       read(*,*)domain !testing the domain 
       if(domain.EQ.1)then
       ox=327962.147  !wf
       oy=2774084.387
       oz=766.
       ox2=270176.646 ! cks  
       oy2=2773555.278
       oz2=27.
        rnorth=80.
        south= 20.
        east=  15.
        west= -75.
       elseif(domain.EQ.2)then
       ox=327962.147  !wf
       oy=2774084.387
       oz=766.
       ox2=313085.626  ! hl
       oy2=2654085.603
       oz2=63.
        west=20.
        east=100.
        south=-120.
        rnorth=-10.
       elseif(domain.EQ.3)then
       ox=299701.532  ! GI 
       oy=2505893.566
       oz=284.        
       ox2=313085.626  ! hl
       oy2=2654085.603
       oz2=63.
        west=25.
        east=85.
        south=30.
        rnorth=100.
       elseif(domain.EQ.4)then
       ox=234181.303  ! kt 
       oy=2422789.221
       oz=41.
       ox2=299701.532  ! GI 
       oy2=2505893.566
       oz2=284.        
        west=60.
        east=130.
        south=-20.
        rnorth=70.
       elseif(domain.EQ.5)then
       ox=156399.173  ! cg
       oy=2560832.205
       oz=38.
       ox2=234181.303  ! kt 
       oy2=2422789.221
       oz2=41.
        west=-50.
        east=10.
        south=-130.
        rnorth=-70.
       elseif(domain.EQ.6)then
       ox=156399.173  ! cg
       oy=2560832.205
       oz=38.
       ox2=110586.785  ! qc 
       oy2=2607290.13
       oz2=48.
        west=-110.
        east=-30.
        south=-60.
        rnorth=20.
       elseif(domain.EQ.7)then
       ox=156399.173  ! cg
       oy=2560832.205
       oz=38.
       ox2=212865.285 ! mq 
       oy2=2682867.833
       oz2=203.
        west=-70.
        east=0.
        south=60.
        rnorth=130.
       elseif(domain.EQ.8)then
       ox=270176.646 ! cks  
       oy=2773555.278
       oz=27.
       ox2=212865.285 ! mq 
       oy2=2682867.833
       oz2=203.
        west=-110.
        east=-40.
        south=-30.
        rnorth=30.
       endif
       nzmax=40 !250m
c       nzmax=20 !500m
       nxmax=int(east-west)+1
       nymax=int(rnorth-south)+1

       print *,'nxmax=',nxmax,'nymax=',nymax 
c      close(10)
       call plot_dual(ixmax,iymax,filein,
     &   ox,oy,oz,ox2,oy2,oz2,south,rnorth,west,east,nxmax,nymax,nzmax,
     &   domain)
       end
       
       subroutine plot_dual(ixmax,iymax,filein,
     &   ox,oy,oz,ox2,oy2,oz2,south,rnorth,west,east,nxmax,nymax,nzmax,
     &   domain)
       real uf(ixmax,iymax),vf(ixmax,iymax),dbz(ixmax,iymax)
       real wfv(ixmax,iymax),vtv(ixmax,iymax)
       real vr_radar1(ixmax,iymax),vr_radar2(ixmax,iymax)
       real dbz_clv(21)
       real wfv_clv(16),wfv1_clv(16),wfv2_clv(16)
       integer pmode
       real*8 gaea_x,gaea_y
       real rx(3000),ry(3000)
       real ty_speed,ty_u,ty_v,ty_theta  !for typhoon rainband
       real px4,py4
       integer domain,x_center,y_center
       integer iasf(13),ity(21)
       character llb(16)*2,ca7*14,ca8*14,ca9*14,date*12
       character ca2*2,ca4*4,ca1*1,ca3*3,ca5*5,case_name*8
       character filein*100,cyy*4,cmo*2,cdd*2,chh*2,cmm*2,clvl*5
       common /pattern/ity
       data iasf/13*1/
       data ity/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21/
       real external Vt  ! terminal velocity calc
       real external variance
       real*8 xdis_rd(2),ydis_rd(2),zdis_rd(2),R(2),vr_vari(2)
       real*8 Sxx,Syy,Sxy,Sxz,Syz,tmp,vtmp,wtmp,kapa
       real*8 w_vari(nxmax,nymax,nzmax+1)
       real*8 u_vari(nxmax,nymax,nzmax)
       real*8 v_vari(nxmax,nymax,nzmax)
       real plot_w(nxmax,nymax),plot_u(nxmax,nymax),plot_v(nxmax,nymax)
       real z_height(nzmax),plot_vv(nzmax)
       real   cont_index(32)
       character dm_ti*10,lv_ti*10
C-------------- pattern center -----------------------
       ca7='$$$$$$$$$$$$$$'
       ca8="$'$'$'$'$'$'$'"
       ca9="$$$$'$'$$$$'$'"
       pi=4.*atan(1.)
       ctr=pi/180.
       dummy=-1000.


c=====calculate vt using height and dBZ
cc     calculate variance of vr,w,vt, define scale height Hs

       vr_vari(1) = 1. ! radial v m/s
       vr_vari(2) = 1. ! radial vm/s
       vt_vari = 1.   ! terminal v m/s
       w_vari = 0. !initial w_vari= 0. m/s
       Hs = 9581.25 !scale height from Gray's inner 2deg composite
       delta_z=500. ! vertical grid spacing (meters)
       delta_x=1000. ! x grid spacing (meters)
       delta_y=1000. ! y grid spacing (meters)
c       kapa=1./Hs
c=====calculate variance==============

      do 310 nxx=int(west),int(east)      
      do 311 nyy=int(south),int(rnorth)
c==== for array index
      ipx=nxx-int(west)+1
      ipy=nyy-int(south)+1
c================================================
cc     the point where the radial wind is observed
      ptx=float(nxx)
      pty=float(nyy)
      do 304 k=nzmax,1,-1  ! nlevel

      pheight=float(k*int(delta_z/2.)) !the point height
c      pheight=float(k*250)
      kapa=-1./(exp(-pheight/Hs))*
     &   ((exp(-(pheight+delta_z)/Hs)-exp(-pheight/Hs))/(delta_z))
c      kapa=-1./(exp(-pheight/Hs))*
c     &    ((exp(-(pheight+500)/Hs)-exp(-(pheight-500)/Hs))/(2*delta_z))
      
c      print *,kapa
      pox=ox + (ptx*1000.)
      poy=oy + (pty*1000.)
c     calc the distance x y z
c     xdis_rd == x-xi , ydis_rd == y-yi , zdir_rd == z-zi
      xdis_rd(1) = (pox-ox)  ! distance btw rd1 and pt (meters)
      ydis_rd(1) = (poy-oy) 
      zdis_rd(1) = (pheight-oz)
      xdis_rd(2) = (pox-ox2) ! distance btw rd2 and pt (meters)
      ydis_rd(2) = (poy-oy2) 
      zdis_rd(2) = (pheight-oz2)

c      print *,xdis_rd,ydis_rd,zdis_rd,'xyz'
c     R == Ri  R(1)=radar1, R(2)=radar2
c     calc Sxx, Syy, Sxy, Sxz, Syz, R
      Sxx = 0.
      Syy = 0.
      Sxy = 0.
      Sxz = 0.
      Syz = 0.
      do i=1,2
       R(i) = (xdis_rd(i)**2+ydis_rd(i)**2+zdis_rd(i)**2)**0.5
       Sxx = Sxx + xdis_rd(i)**2
       Syy = Syy + ydis_rd(i)**2
       Sxy = Sxy + xdis_rd(i)*ydis_rd(i)
       Sxz = Sxz + xdis_rd(i)*zdis_rd(i)
       Syz = Syz + ydis_rd(i)*zdis_rd(i)
      enddo
c=====calculate variance of u,v,w
      do 307 kk=1,100
       utmp=0.
       vtmp=0.
c     calculate u_variance , v_variance
       do i=1,2
        utmp = utmp +
     &    (vr_vari(i)*(R(i)*xdis_rd(i)*Syy-R(i)*ydis_rd(i)*Sxy)**2)
     &    /((Sxx*Syy-(Sxy)**2)**2)
        vtmp = vtmp +
     &    (vr_vari(i)*(R(i)*ydis_rd(i)*Sxx-R(i)*xdis_rd(i)*Sxy)**2)
     &    /((Sxx*Syy-(Sxy)**2)**2)
c       print *,utmp,vtmp,'utmp vtmp'
       enddo
c      variance u at k level
c      variance v at k level
       u_vari(ipx,ipy,k) = utmp+
     &         ((w_vari(ipx,ipy,k)+vt_vari)*((Sxy*Syz-Syy*Sxz))**2
     &        /((Sxx*Syy-(Sxy)**2))**2)
       v_vari(ipx,ipy,k) = vtmp+
     &         ((w_vari(ipx,ipy,k)+vt_vari)*((Sxy*Sxz-Sxx*Syz))**2
     &        /((Sxx*Syy-(Sxy)**2))**2)

c      calculate w_variance
       w_vari(ipx,ipy,k)=
     &         0.5*(1+((1/delta_z-kapa/2.)/(1/delta_z+kapa/2.))**2)
     &       * w_vari(ipx,ipy,k+1)
     &       + u_vari(ipx,ipy,k)/((2.*delta_x*(1./delta_z+kapa/2.))**2)
     &       + v_vari(ipx,ipy,k)/((2.*delta_y*(1./delta_z+kapa/2.))**2)
307   enddo
cc     write down the result 
c      write(15,*)'level',pheight,
c     &         'u_vari',u_vari(k),'v_vari',v_vari(k),'w_vari',w_vari(k)
c=====end of codes====================
304   enddo

311   enddo
310   enddo

c      do i=1,20
c       write(15,*)'level',i*500,
c     &    u_vari(46,41,i),v_vari(46,41,i),w_vari(46,41,i)
c      enddo
      
      call opngks
      call defclr16
      call setusv('LW',6000)
      call gsclip(1)

      if(nymax.gt.nxmax)then
       dd=0.8*(float(nxmax-1)/float(nymax-1))

       call set((0.5-dd/2.),(0.5+dd/2.),0.1,0.9,
     +           1.,float(nxmax),1.,float(nymax),1)

       elseif(nymax.eq.nxmax)then
       dd=0.8*(float(nxmax-1)/float(nymax-1))
       call set(0.1,0.1+dd,0.1,0.9,
     +           1.,float(nxmax),1.,float(nymax),1)
       else
       dd=0.8*(float(nymax-1)/float(nxmax-1))
       call set(0.1,0.9,0.1,0.1+dd,
     +           1.,float(nxmax),1.,float(nymax),1)
       endif
       call getset(f1,f2,f3,f4,u1,u2,u3,u4,l1)
c-----------------------------------------
       write(*,*)'input LEVEL(m)='
       read(*,*)nlevel
c       nlevel=1000
       nlevel=nlevel/int(delta_z/2)
       print *,'nlevel',nlevel
       do i=1,nxmax
        do j=1,nymax
         plot_w(i,j)=w_vari(i,j,nlevel)
        enddo
       enddo
c      print *,'~~~~~'
       where(plot_w.gt.100.)plot_w=dummy
c       do i=1,nxmax
c         write(15,*)(plot_w(i,j),j=1,nymax)
c       enddo

       do i=1,32
        cont_index(i)=0.0+float(i-1)*0.2
       enddo
      call gsplci(16)
c-------plot vertical velocity error---------------
      call contour(plot_w,nxmax,nymax,dummy,cont_index,32,1.5,ca7,16,4)
      call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      call set(f1,f2,f3,f4,u1,u2,u3,u4,l1)
      idx=int((east-west)/10.)
      idy=int((rnorth-south)/10.)
      call axis(idx,1,'X (km)',idy,1,'Y (km)')
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(dm_ti,'(i2)')int(domain)
      call plchhq(0.5,ft+0.04,'w error Domain '//trim(dm_ti),0.02,0.,0.)
      write(lv_ti,'(f4.1)')float(nlevel)*delta_z/2000.
      call plchhq(0.8,ft+0.04,'level'//trim(lv_ti)//'km',0.014,0.,0.)
      call frame

c------plot u component error-----------------------
       do i=1,nxmax
        do j=1,nymax
         plot_u(i,j)=u_vari(i,j,nlevel)
        enddo
       enddo
      call set(f1,f2,f3,f4,u1,u2,u3,u4,l1)
      call contour(plot_u,nxmax,nymax,dummy,cont_index,32,1.5,ca7,16,4)
      call gsplci(16)
      call set(fl,fr,fb,ft,west,east,south,rnorth,l1)
      idx=int((east-west)/10.)
      idy=int((rnorth-south)/10.)
      call axis(idx,1,'X (km)',idy,1,'Y (km)')
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(dm_ti,'(i2)')int(domain)
      call plchhq(0.5,ft+0.04,'u error Domain '//trim(dm_ti),0.02,0.,0.)
      write(lv_ti,'(f4.1)')float(nlevel)*delta_z/2000.
      call plchhq(0.8,ft+0.04,'level'//trim(lv_ti)//'km',0.014,0.,0.)
      call frame

c------plot v component error-----------------------
       do i=1,nxmax
        do j=1,nymax
         plot_v(i,j)=v_vari(i,j,nlevel)
        enddo
       enddo
      call set(f1,f2,f3,f4,u1,u2,u3,u4,l1)
      call contour(plot_v,nxmax,nymax,dummy,cont_index,32,1.5,ca7,16,4)
      call gsplci(16)
      call set(fl,fr,fb,ft,west,east,south,rnorth,ll)
      idx=int((east-west)/10.)
      idy=int((rnorth-south)/10.)
      call axis(idx,1,'X (km)',idy,1,'Y (km)')
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(dm_ti,'(i2)')int(domain)
      call plchhq(0.5,ft+0.04,'v error Domain '//trim(dm_ti),0.02,0.,0.)
      write(lv_ti,'(f4.1)')float(nlevel)*delta_z/2000.
c      print *,lv_ti
      call plchhq(0.8,ft+0.04,'level'//trim(lv_ti)//'km',0.014,0.,0.)
      call frame

c--------plot value with height in the center of domain--------
      x_center=int((west+east)/2.)-west+1
      y_center=int((south+rnorth)/2.)-south+1
      print *,x_center,y_center
      do i=1,3
       plot_vv=dummy
       do k=1,nzmax
        z_height(k)=250.+250*float(k-1)
        select case(i)
         case(1)
         plot_vv(k)=w_vari(x_center,y_center,k)
         case(2)
         plot_vv(k)=v_vari(x_center,y_center,k)
         call dashdb(130175)
         case(3)
         plot_vv(k)=u_vari(x_center,y_center,k)
         call dashdb(52428)
        end select
        write(21,*)z_height(k),plot_vv(k)
       enddo
      call set(0.1,0.9,0.1,0.9,0.,4.,1000.,10000.,1)
      call curved(plot_vv,z_height,nzmax)
      enddo
      call axis(4,1,'variance (m:S:-2:N:s:S:-2:N:)',18,1,'Height (km)')
      call frame
      call clsgks

      END



      subroutine lbar(llb,n1,n2,mode)
      dimension iasf(13)
      integer ll, n1, n2
      real fl,fr,ft,fb,ul,ur,ut,ub
      character*2 llb(16)
      integer lnd(16)
      common /pattern/ lnd
      data iasf/13*1/
        call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
        call gsclip(1)
        call gsasf(iasf)
c        call gsfais(1)

c       call plotif(0.,0.,2)
       call gsplci(16)
       call gstxci(16)
       call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c       call lblbar(0,0.1,0.9,.92,0.99,n1,1.,0.4,lnd,0,llb,n2,1)
       call lblbar(0,fl,fr,ft+0.02,ft+0.09,n1,1.,0.4,lnd,0,llb,n2,1)

       if(mode .eq. 1)then
       call plchhq(0.5,ft+0.015,'Vertical velocity (m/s) ',
     +  0.015,0.,0.)
       elseif(mode .eq. 2)then
        call plchhq(0.5,ft+0.015,'Radar Reflectivity (dBZ) ',
     +  0.015,0.,0.)
       endif
       call set(fl,fr,fb,ft,ul,ur,ub,ut,ll)
       return
       end


      SUBROUTINE color(dat,nx,ny,clv,dummy)
       integer nx,ny
       dimension dat(nx,ny)
c       dimension iama(150000),rwrk(70000),iwrk(40000)    !by yuku
       dimension iama(1600000),rwrk(1600000),iwrk(1600000) ! by lwcheng
c       dimension xcra(8500),ycra(8500),iara(10),igra(10) !by yuku
       dimension xcra(160000),ycra(160000),iara(10),igra(10) !by lwcheng
       dimension clv(21)
       dimension iasf(13)

       EXTERNAL colram
       data iasf /13*1/

       call gsclip(1)
       call gsasf(iasf)
       call gsfais(1)
       call SFSETI('TYPE OF FILL',0)
       call SFSETI ('ANGLE OF FILL LINE',0)
       call SFSETR('SPACING OF FILL LINES',.001)

       call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
       call CPSETR('VPS',0.)
       call CPSETR('VPL - VIEWPORT LEFT',fl) 
       call CPSETR('VPR - VIEWPORT RIGHT',fr)
       call CPSETR('VPB - VIEWPORT BOTTOM',fb) 
       call CPSETR('VPT - VIEWPORT TOP',ft) 
       call CPSETR('SPV - SPECIAL VALUE',dummy)
cc       call CPSETR('T2D - TENSION ON 2-DIMENSION SPLINES',4.)
       call CPSETI('NOF - NUMERIC OMISSION FLAGS', 7)
       call CPSETI('CLS - CONTOUR LEVEL SELECTION', 0)
       call CPSETI('NCL - NUMBER OF CONTOUR LEVELS',21)
       DO I = 1,20
          call CPSETI('PAI - PARAMETER ARRAY IDENTIFIER',I)
       call CPSETR('CLV - CONTOUR LEVEL VALUE',clv(I))
       call CPSETI('CLU - CONTOUR LEVEL USE',1)
c	  call CPSETR('CLL - CONTOUR LEVEL LINE WIDTH',2.)
       call CPSETI('AIA - AREA IDENTIFIER ABOVE LINE',I+1)
       call CPSETI('AIB - AREA IDENTIFIER BELOW LINE',I  )
       END DO 
cc       call CPRECT(dat,nx,nx,ny,RWRK,70000,IWRK,40000)  !by yuku
       call CPRECT(dat,nx,nx,ny,RWRK,1600000,IWRK,1600000)  !by lwcheng
c       call ARINAM(IAMA,150000) !by yuku
       call ARINAM(IAMA,1600000) !by lwcheng
       call CPCLAM(dat,RWRK,IWRK,IAMA)  ! error
       call cpgeti('IWU',iiwu)
       call cpgeti('RWU',irwu)
       call cpgeti('WSO',iwso)
       write(*,*) 'IWU RWU WSO ',iiwu,irwu,iwso
c       call ARSCAM(IAMA,XCRA,YCRA,8500,IARA,IGRA,10,colram) !by yuku
       call ARSCAM(IAMA,XCRA,YCRA,1600000,IARA,IGRA,10,colram) !by lwcheng
c       call CPSETC('ILT','')
c       call CPLBDR(dat,RWRK,IWRK)
c       call CPCLDR(dat,RWRK,IWRK)

       call SET(fl,fr,fb,ft,ul,ur,ub,ut,ll)
       RETURN
       END

C----------------------------------------------------------------------
C     The arrays XCRA and YCRA, for indices 1 to NCRA, contain the X
C     and Y coordinates of points defining a polygon. The area
C     identifier in the array IAIA, each with an associated group
C     identifier in the array IGIA, tell us whether the polygon is to 
C     be color-filled or not.
C----------------------------------------------------------------------
      SUBROUTINE COLRAM(XCRA,YCRA,NCRA,IAIA,IGIA,NAIA)

      DIMENSION XCRA(*),YCRA(*),IAIA(*),IGIA(*),
     +          DST(5000),IND(5000)
       DIMENSION ity(21)
       COMMON /pattern/ ity

C Assume the polygon will be filled until we find otherwise.
      NST = 5000
      NND = 5000
      IFLL = 1

C If any of the area identifiers is nagative, don't fill the polygon.

      DO 101 I = 1, NAIA
         IF (IAIA(I) .LT. 0 ) IFLL = 0 
  101 CONTINUE

C Otherwise, fill the polygon in the color implied by its area
C identifier relative to edge group 3 (the contour-line group).
 
      IF (IFLL .NE. 0 ) THEN
         IFLL = 0
         DO I = 1, NAIA
           IF (IGIA(I) .EQ. 3 ) IFLL = ity(IAIA(I)-1)
         END DO
         IF ( IFLL .GT. 0 .AND. IFLL .LT. 20 ) THEN
C	    CALL GSFACI(IFLL+1)
C	    CALL GFA( NCRA, XCRA, YCRA)
C	    CALL SFSETI('TY',1)
            CALL SFSGFA(XCRA, YCRA, NCRA, DST,NST, IND,NND, IFLL)
         END IF
      END IF

      RETURN
      END

      SUBROUTINE DEFCLR16
      DIMENSION RGBV(3,15)
      DATA RGBV /
     +            0.70,0.70,0.70,
     +            0.00,1.00,1.00,
     +            0.00,0.75,1.0,
     +            0.00,0.50,1.0,
     +            0.00,0.50,0.30,
     +            0.00,0.75,0.30,
     +            0.00,1.00,0.30,
     +            1.00,1.00,0.00,
     +            1.00,0.75,0.00,
     +            1.00,0.00,0.00,
     +            0.50,0.00,0.00,
     +            0.81,0.34,0.64,
     +            1.00,0.50,0.00,
     +            0.40,0.40,0.40,
     +            0.85,0.63,0.80/
c      DATA RGBV / 0.90,0.90,0.90,
c     +            0.60,1.00,0.00,
c     +            0.00,0.85,0.90,
c     +            0.00,0.60,0.80,
c     +            0.00,0.30,0.70,
c     +            0.50,0.00,0.00,
c     +            0.35,0.10,0.10,
c     +            0.10,0.20,0.10,
c     +            0.00,0.75,0.00,
c     +            1.00,0.38,0.38,
c     +            1.00,0.00,0.38,
c     +            1.00,0.00,0.00,
c     +            0.00,0.00,0.00,
c     +            0.00,0.00,0.00,
c     +            1.00,1.00,1.00 /

c       do i=1,9
c       dd=0.11*float(i)
c       RGBV(1,i)=1.-dd
c       RGBV(2,i)=1.-dd
c       RGBV(3,i)=1.-dd
c       enddo

       CALL GSCR (1,0,1.,1.,1.)  !    background color (1,0,0.,0.,0.)<--black 
       DO 101 I = 1, 15
       CALL GSCR (1,I,RGBV(1,I),RGBV(2,I),RGBV(3,I))
  101  CONTINUE

       CALL GSCR (1,16,0.,0.,0.)
       CALL GSCR (1,22,0.55,0.44,0.00)
       RETURN
       END



      SUBROUTINE FLAG(X,Y,SPEED,THETA,IUV,LEN)
C
C IUV=1 (SPEED,THETA) PLOT FLAG
C    =2 (U,V) PLOT FLAG
C    =-1 (SPEED,THETA) PLOT VECTOR
C    =-2 (U,V) PLOT VECTOR
C
      PARAMETER (SWP=7.5/11., DXL=0.2, SPDMIN=3.)
      PARAMETER (DEG120=120.,N2=2, PI=3.1415926)
      PARAMETER (V10=10., HFLN=0.5)
      PARAMETER (N3=N2+1, DEG60=PI-PI/180.*DEG120)
      DATA SPV/-999./
      CALL GETSET(r11,r12,r13,r14,X11,X22,Y11,Y22,N15)
      n11=1024*r11
      n12=1024*r12
      n13=1024*r13
      n14=1024*r14
      RLENG=ABS(LEN)
      IF(LEN.EQ.0)RLENG=25.
      RLENGX=RLENG/(N12-N11)*(X22-X11)
      RLENGY=RLENG/(N14-N13)*(Y22-Y11)
      IF(IABS(IUV) .LE. 1)THEN
          THE=(90. - THETA) * PI / 180.
          SPD=SPEED
          IF(SPD.LE.0.)GOTO 50
      ELSE
          U=SPEED
          V=THETA
          SPD=SQRT( U*U + V*V )
          IF(SPD .EQ. 0.)GOTO 50
          THE=ATAN2( -V, -U)
      ENDIF
      IF(SPD.EQ.SPV)GOTO 50
      IF(LEN.LT.0)SPD=SPD*2
      IF(SPD.EQ.SPV)GOTO 50
      CALL FRSTPT(X,Y)
      IF(IUV.LT.0)GOTO 30
C
C  --- PLOT WIND FLAG
C
      SP120X=RLENGX * SWP
      SP120Y=RLENGY * SWP
      XL= RLENGX * COS(THE)
      YL= RLENGY * SIN(THE)
      DLX = XL * DXL
      DLY = YL * DXL
      IF(SPD.GE.83)THEN
        XL=XL*1.5
        YL=YL*1.5
      ENDIF
      X2=X+XL
      Y2=Y+YL
      CALL VECTOR(X2,Y2)
      IF(SPD .LT.SPDMIN)GOTO 50
      D1=THE - DEG60
      XE= SP120X * COS(D1)
      YE= SP120Y * SIN(D1)
              IF(SPD.GE.8.)THEN
C
C  --- 50 KNOT
C
      LS1= (SPD+2) / 50
      DO 10 I=1,LS1
      X22=X2-DLX+XE
      Y22=Y2-DLY+YE
      CALL VECTOR(X22,Y22)
      DO 8 K=1,N2
      AX=DLX*K/N3
      AY=DLY*K/N3
      CALL FRSTPT(X2-AX,Y2-AY)
   8  CALL VECTOR(X22,Y22)
      X2=X2 - DLX
      Y2=Y2 - DLY
      CALL VECTOR(X2,Y2)
      X2=X2 - DLX
      Y2=Y2 - DLY
      CALL FRSTPT(X2,Y2)
  10  CONTINUE
      REM=SPD-LS1*50
      IF(REM .LT. 3.)GOTO 50
C
C  --- 10 KNOT
C
      LS2= (REM+2) / 10
      DO 20 I=1,LS2
      CALL VECTOR(X2+XE,Y2+YE)
      X2=X2 - DLX
      Y2=Y2 - DLY
      CALL FRSTPT(X2,Y2)
  20  CONTINUE
C
C  --- 5 KNOT
C
      REM= SPD - INT((SPD+2)/10) * 10
      IF(REM .LT. 3.)GOTO 50
              ELSE
      X2=X + XL * 0.6
      Y2=Y + YL * 0.6
      CALL FRSTPT(X2,Y2)
              ENDIF
      CALL VECTOR(X2+XE/2, Y2+YE/2)
      GOTO 50
  30  CONTINUE
C
C  --- PLOT WIND VECTOR
C
      THE=THE+PI
      XY=SPD/V10
      X2=X+XY*COS(THE)*RLENGX
      Y2=Y+XY*SIN(THE)*RLENGY
      CALL VECTOR(X2,Y2)
      THE1=THE+PI/6
      X1=-COS(THE1)*RLENGX*HFLN*XY
      Y1=-SIN(THE1)*RLENGY*HFLN*XY
      CALL VECTOR(X2+X1,Y2+Y1)
      THE1=THE-PI/6
      X1=-COS(THE1)*RLENGX*HFLN*XY
      Y1=-SIN(THE1)*RLENGY*HFLN*XY
      CALL FRSTPT(X2,Y2)
      CALL VECTOR(X2+X1,Y2+Y1)
  50  CONTINUE
C      CALL FRSTPT(X,Y)
      RETURN
      END








      subroutine get_topo_1000(ox,oy,west,south,topo,ix,iy,mode)
C-----------------------------------------------------------------
C ------------   mode: useing interpolation or not (1:no , 2:yes)
C    writen by lwcheng(2005.9.11)
C-----------------------------------------------------------------
      real    topo(ix,iy)
      integer hei(201,377)

c       gaea_x=150555.7869           
c       gaea_y=2422360.3752         
       rox=(ox - 150555.-0.7869)/1000.
       roy=(oy- 2422360.-0.3752)/1000.
       rox=rox+west
       roy=roy+south


       open(20,file='/home/lwcheng/DTM/linux_X1000m.hei'
     +              ,form='unformatted')
       read(20)hei
       close(20)

       iix=rox
       iiy=roy

       if(mode.eq.1)then
       do j=1,iy
       do i=1,ix
       topo(i,j)=0.
         ii=iix+i
         jj=iiy+j

       if(ii.ge.1.and.ii.le.201.and.jj.ge.1.and.jj.le.377)
     + topo(i,j)=hei(ii,jj)

       enddo
       enddo


C-----------------------  interpolation of topo  -----------------
       elseif(mode.eq.2)then
       do j=1,iy
       do i=1,ix
       topo(i,j)=0.
       rrox=rox+float(i)
       rroy=roy+float(j)
       iix=rrox
       iiy=rroy
       ii=iix
       jj=iiy


c       if(ii.ge.1.and.ii.le.200.and.jj.ge.1.and.jj.le.376)then
       if(ii.ge.1.and.ii.le.200.and.jj.ge.1.and.jj.le.376
     + .and. hei(ii,jj)   .ne.0 .and. hei(ii+1,jj)  .ne.0 .and.
     +       hei(ii,jj+1) .ne.0 .and. hei(ii+1,jj+1).ne.0 )then

       v1=hei(ii,jj)+(hei(ii+1,jj)-hei(ii,jj))*(rrox-float(iix))
       v2=hei(ii,jj+1)+(hei(ii+1,jj+1)-hei(ii,jj+1))*(rrox-float(iix))
       topo(i,j)=v1+(v2-v1)*(rroy-float(iiy))
       endif

       enddo
       enddo

       endif

      return
      end


      subroutine get_topo_400(ox,oy,west,south,topo,ix,iy,mode)
C------------------------------------------------------------------------ 
C ------------   mode: useing interpolation or not (1:no , 2:yes)
C        [1000m]  => [400m]
C        formula  => (i-1)*2.5+1
C   topo(181,181) => (i-1)*2.5+1 =>  topo(451,451)
C   writen by lwcheng(2005.9.11)
C-------------------------------------------------------------------------
      real topo(ix,iy)
      integer hei(502,943)

c       gaea_x=150555.7869           
c       gaea_y=2422360.3752         
       rox=(ox - 150555. -0.7869 +(west*1000.) )/400.
       roy=(oy- 2422360. -0.3752 +(south*1000.))/400.
c       ox=ox+(west)40.
c       oy=oy+(south)40.


       open(20,file='/home/lwcheng/DTM/linux_X400m.hei'
     +        ,form='unformatted')
       read(20)hei
       close(20)

       iix=rox
       iiy=roy

       if(mode.eq.1)then
       do j=1,iy
       do i=1,ix
       topo(i,j)=0.
         ii=iix+i
         jj=iiy+j

       if(ii.ge.1.and.ii.le.502.and.jj.ge.1.and.jj.le.943)
     + topo(i,j)=hei(ii,jj)

       enddo
       enddo


C-----------------------  interpolation of topo  -----------------
       elseif(mode.eq.2)then
       do j=1,iy
       do i=1,ix
       topo(i,j)=0.
       rrox=rox+float(i)
       rroy=roy+float(j)
       iix=rrox
       iiy=rroy
       ii=iix
       jj=iiy


c       if(ii.ge.1.and.ii.le.501.and.jj.ge.1.and.jj.le.942)then
       if(ii.ge.1.and.ii.le.501.and.jj.ge.1.and.jj.le.942
     + .and. hei(ii,jj).ne.0 .and. hei(ii+1,jj).ne.0 .and.
     +       hei(ii,jj+1).ne.0 .and. hei(ii+1,jj+1).ne.0 )then

       v1=hei(ii,jj)+(hei(ii+1,jj)-hei(ii,jj))*(rrox-float(iix))
       v2=hei(ii,jj+1)+(hei(ii+1,jj+1)-hei(ii,jj+1))*(rrox-float(iix))
       topo(i,j)=v1+(v2-v1)*(rroy-float(iiy))
       endif

       enddo
       enddo

       endif

      return
      end

      subroutine get_topo_40(ox,oy,west,south,topo,ix,iy)
C-------------------------------------------------------------
C        [1000m]   =>   [40m]
C        formula   =>   (i-1)*25+1
C    topo(51,51)   =>   (i-1)*25+1 =>  topo(1251,1251) 
C    writen by lwcheng(2005.9.11)
C-------------------------------------------------------------
      real topo(ix,iy)
      integer hei(5012,9425)
c       gaea_x=150555.7869           
c       gaea_y=2422360.3752         
       rox=(ox - 150555. -0.7869 +(west*1000.) )/40.
       roy=(oy- 2422360. -0.3752 +(south*1000.))/40.


       open(20,file='/home/lwcheng/DTM/linux_X40m.hei'
     +        ,form='unformatted')
       read(20)hei

       iix=rox
       iiy=roy

       do j=1,iy
         do i=1,ix
           topo(i,j)=0.
           ii=iix+i
           jj=iiy+j
             if(ii.ge.1.and.ii.le.5012.and.jj.ge.1.and.jj.le.9425)then
             topo(i,j)=hei(ii,jj)
             endif
         enddo
       enddo
       close(20)

      return
      end

      subroutine shader_topo(ox,oy,east,west,south,rnorth,iclr)
      real topo(10001,10001),clv_hei(8)
      real rx(4000),ry(4000)
      real*8 gaea_x,gaea_y
      character llb(16)*4,ca7*14,ca8*14,ca9*14
      integer lsha(16)
      integer iasf(13)
      data iasf/13*1/
      data lsha/0,6,8,10,12,14,14,14,14,10,11,12,13,14,15,16/
c      data lsha/0,6,8,10,5,6,7,8,9,10,11,12,13,14,15,16/
      common /shad/lsha
      ca7='$$$$$$$$$$$$$$'
      ca8="$'$'$'$'$'$'$'"
      ca9="$$$$'$'$$$$'$'"
      call gsplci(iclr)
 
      ix=nint((east-west)*2.5)+1
      iy=nint((rnorth-south)*2.5)+1
      call get_topo_400(ox,oy,west,south,topo(1:ix,1:iy),ix,iy,2)

      
      call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
   
      clv_hei(1)=100.
      clv_hei(2)=500.
      clv_hei(3)=1500.
      clv_hei(4)=2500.
      clv_hei(5)=4000.
      clv_hei(6)=5000.
      clv_hei(7)=6000.
      clv_hei(8)=7000.


      do i=1,6
      write(llb(i),'(I4)')int(clv_hei(i))
      enddo
      llb(6)='    '


C------------------ topo shader --
      if(1.eq.2)then

      call setusv("LW",2000)
      call shader(topo(1:ix,1:iy),ix,iy,clv_hei(2:6),5,-999.)

      lsha(1)=6
      call lbar2(llb,3)
      lsha(1)=8
      call lbar2(llb,4)
      lsha(1)=10
      call lbar2(llb,5)


      call lbar3(llb)

      call contour(topo(1:ix,1:iy),ix,iy,0.,clv_hei(1),1,1.,ca7,iclr,3)
      call contour(topo(1:ix,1:iy),ix,iy,0.,clv_hei(2),1,1.,ca7,iclr,5)
      call contour(topo(1:ix,1:iy),ix,iy,0.,clv_hei(3),1,1.,ca7,iclr,3)
      call contour(topo(1:ix,1:iy),ix,iy,0.,clv_hei(4),1,1.,ca7,iclr,3)
      call contour(topo(1:ix,1:iy),ix,iy,0.,clv_hei(5),1,1.,ca7,iclr,3)
     
      endif

C---------------- plot tw outline ----------------------
      call setusv('LW',6000) 
      if(ix.gt.iy)then
      dd=float(iy-1)/float(ix-1)*0.8
      call set(0.1,0.9,0.1,0.1+dd,0.,float(ix-1),0.,float(iy-1),1)
      else
      dd=float(ix-1)/float(iy-1)*0.8
      call set(0.1,0.1+dd,0.1,0.9,0.,float(ix-1),0.,float(iy-1),1)
      endif


c       open(33,file='/home/lwcheng/DTM/Xoutline.dat',status='old')
c       open(33,file='/home/lwcheng/DTM/Xtw_coast_city.dat',status='old')
      call gsclip(1)
       open(33,file='/home/lwcheng/DTM/Xtw_coast.dat',status='old')
197   continue
      nn=0
      read(33,*,end=198)number
      do i=1,number
      read(33,*,end=198)gaea_x,gaea_y
      gaea_x=(gaea_x-((ox+west*1000.) ))/400.  
      gaea_y=(gaea_y-((oy+south*1000.) ))/400.
      xx=gaea_x
      yy=gaea_y
      nn=nn+1
      rx(nn)=xx
      ry(nn)=yy
      enddo
c      call curve(rx,ry,nn)
      goto 197
198   continue


      call SET(fl,fr,fb,ft,ul,ur,ub,ut,ll)


      return
      end


      SUBROUTINE shader(dat,nx,ny,clv,n,dummy)
       dimension dat(nx,ny),dat2(nx,ny)
       dimension iama(15000000),rwrk(8000000),iwrk(8000000)
       dimension xcra(7000000),ycra(7000000),iara(10),igra(10)
       dimension clv(n)
       dimension iasf(13)
       data iasf /13*1/
       EXTERNAL SHADAM
       dat2=dat
       where(dat2.lt.clv(1))dat2=dummy
       where(dat2.gt.clv(n))dat2=dummy

       call gsclip(1)
       call gsasf(iasf)
       call SFSETI('TYPE OF FILL',-2)
       call SFSETI ('ANGLE OF FILL LINE',45)
       call SFSETR('SPACING OF FILL LINES',.001)

       call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
       call CPSETR('VPL - VIEWPORT LEFT',fl) 
       call CPSETR('VPR - VIEWPORT RIGHT',fr)
       call CPSETR('VPB - VIEWPORT BOTTOM',fb)
       call CPSETR('VPT - VIEWPORT TOP',ft) 
       call CPSETR('SPV - SPECIAL VALUE',dummy)
C       call CPSETR('T2D - TENSION ON 2-DIMENSION SPLINES',1.)
       call CPSETI('NOF - NUMERIC OMISSION FLAGS', 7)
       call CPSETI('CLS - CONTOUR LEVEL SELECTION', 0)
       call CPSETI('NCL - NUMBER OF CONTOUR LEVELS',n)
       DO I = 1,n-1
          call CPSETI('PAI - PARAMETER ARRAY IDENTIFIER',I)
          call CPSETR('CLV - CONTOUR LEVEL VALUE',clv(I))
c         call CPSETI('CLU - CONTOUR LEVEL USE',3) ! contour and word
          call CPSETI('CLU - CONTOUR LEVEL USE',0) ! no
c         call CPSETR('CLL - CONTOUR LEVEL LINE WIDTH',2.)
          call CPSETI('AIA - AREA IDENTIFIER ABOVE LINE',I+1)
          call CPSETI('AIB - AREA IDENTIFIER BELOW LINE',I  )
       END DO 
       call CPRECT(dat2,nx,nx,ny,RWRK,800000,IWRK,800000)
       call ARINAM(IAMA,1500000)
       call CPCLDR(dat2,RWRK,IWRK)
       call CPCLAM(dat2,RWRK,IWRK,IAMA)
       call cpgeti('IWU',iiwu)
       call cpgeti('RWU',irwu)
       call cpgeti('WSO',iwso)
       write(*,*) 'IWU RWU WSO ',iiwu,irwu,iwso
       call ARSCAM(IAMA,XCRA,YCRA,700000,IARA,IGRA,10,SHADAM)

       call SET(fl,fr,fb,ft,ul,ur,ub,ut,ll)
       RETURN
       END

       SUBROUTINE SHADAM(XCS,YCS,NCS,IAI,IAG,NAI) 
C      This version of SHADAM shades the area-map polygon whose edge is
C      defined by the points ((XCS(I),YCS(I)),I=1,NCS) if and only,
C      relative to edge group 3, its area identifier is between 1 and
C      20. The density of the shading increase with the value of the
C      area identifier.

       DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
       DIMENSION DST(500000),IND(500000)
       DIMENSION ity(16)
       COMMON /shad/ ity
C------Turn off shading ---------------------------------------------
       ISH = 0 

C------If the area identifier fo r group 3 is in the right range,
C------Turn on shading.
       DO I = 1, NAI 
          IF (IAG(I) .EQ. 3 .AND. IAI(I) .GE. 1 .AND. IAI(I) .LE. 16)
     +    ISH = ity(IAI(I)-1)
       END DO

C------If shading is turn on, shade the area. The last point of the
C------edge is redundant and may be omitted.
c   ori     IF (ISH.NE.0) CALL SFSGFA(XCS,YCS,NCS,DST,50000,IND,50000,ISH)
       IF (ISH.NE.0) CALL SFSGFA(XCS,YCS,NCS,DST,500000,IND,500000,ISH)   ! lwcheng

       RETURN
       END

        subroutine lbar3(llb)
          dimension iasf(13)
           integer ll, n1, n2
          real fl,fr,ft,fb,ul,ur,ut,ub
        character*4 llb(16)
             integer lnd(16),lsha(16)
            common /shad/ lnd

          lnd=0
          data iasf/13*1/
        call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
        call gsclip(1)
            call gsasf(iasf)
          call gsfais(1)
        call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c       call lblbar(1,0.90,1.,0.1,0.5,n1,0.4,1.,lnd,0,llb,n2,1)
ccccc       call lblbar(1,0.90,0.96,0.4,0.6,2,0.6,1.,lnd,0,llb,3,0)
c        call lblbar(0,0.1,0.9,.75,0.8,n1,1.,0.4,lnd,0,llb,n2,1)

       call plotif(0.,0.,2)
c       call gsplci(16)
c       call gstxci(16)

       call setusv("LW",3000)
       call line(0.913,0.1,0.948,0.1)
       call setusv("LW",5000)
       call line(0.913,0.18,0.948,0.18)
       call setusv("LW",3000)
       call plchhq(0.905,0.555,'Terrain',0.01,0.,-1.)
       call plchhq(0.905,0.535,'Height(m)',0.01,0.,-1.)
       dy=(0.6-0.2)/float(5)
       do i=1,6
       yy=0.1+float(i-1)*dy
       call plchhq(0.975,yy,llb(i),0.01,0.,0.)
       enddo

       call set(fl,fr,fb,ft,ul,ur,ub,ut,ll)
        return
       end


      SUBROUTINE contour(dat,nx,ny,dummy,clv,kn,cwm,cadash,ind,lw)
c**  kn      : contour line no.
c**  cadash  : contour line dash pattern (character*16)
c**  cwm     : character size (real) (default=1.)
c**  ind     : contour line color index (default=16)
c**  lw      : contour line width (default=1) ! by lwcheng 

       integer nx,ny,kn
       dimension dat(nx,ny)
c       dimension rwrk(40000),iwrk(30000) !by yuku
       dimension rwrk(160000),iwrk(160000) !by lwcheng
C       dimension xcra(2500),ycra(2500),iara(10),igra(10)
c       character*16 cadash
       character*14 cadash
       dimension clv(kn)
       dimension iasf(13)

c       EXTERNAL SHADAM
       data iasf /13*1/


       call gsclip(1)
       call gsasf(iasf)

       call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
       call cpseti('LLP -line label',1)
c         cwm = 1. / ( fr - fl )
       call pcseti('QU - quality flag',0)
       call CPSETR('CWM - CHARACTER width Multiplier',cwm)
       call CPSETR('VPL - VIEWPORT LEFT',fl)
       call CPSETR('VPR - VIEWPORT RIGHT',fr)
       call CPSETR('VPB - VIEWPORT BOTTOM',fb)
       call CPSETR('VPT - VIEWPORT TOP',ft)
       call CPSETR('SPV - SPECIAL VALUE',dummy)
       call CPSETR('PC4 - Penalty Scheme Constant 4',.1)
       call CPSETR('PC5 - Penalty Scheme Constant 5',.2)
c       call CPSETR('T2D - TENSION ON 2-DIMENSION SPLINES',4.)
       call CPSETI('NOF - NUMERIC OMISSION FLAGS', 7)
       call CPSETI('NSD - NUMBER of SIGNIFICANT DIGITS', 4)
       call CPSETI('CLS - CONTOUR LEVEL SELECTION', 0)
       call CPSETI('NCL - NUMBER OF CONTOUR LEVELS',kn)
c       call CPSETR('CFS',0.05)
c       call CPSETR('LLS - LINE LABEL SIZE',0.05)
c       call CPSETR('CWM',1.2)
c       call CPSETR('LLW',0.004)
       DO I = 1, kn
          call CPSETI('PAI - PARAMETER ARRAY IDENTIFIER',I)
          call CPSETR('CLV - CONTOUR LEVEL VALUE',clv(I))
          call cpsetc('CLD - CONTOUR LINE DASH PATTERN',cadash)
          call cpseti('CLC - contour line color index',ind)
c         call cpseti('CLC - contour line color index',0)!!!lwcheng
          call CPSETI('CLU - CONTOUR LEVEL USE',3)
c         call CPSETI('CLU - CONTOUR LEVEL USE',1)
          call CPSETR('CLL - CONTOUR LEVEL LINE WIDTH',float(lw))
          call CPSETI('AIA - AREA IDENTIFIER ABOVE LINE',I+1)
          call CPSETI('AIB - AREA IDENTIFIER BELOW LINE',I  )
       END DO
       call CPRECT(dat,nx,nx,ny,RWRK,40000,IWRK,30000) !by yuku
c       call CPRECT(dat,nx,nx,ny,RWRK,160000,IWRK,160000) !by lwcheng
       call CPCLDR(dat,RWRK,IWRK)
       call cpgeti('IWU',iiwu)
       call cpgeti('RWU',irwu)
       call cpgeti('WSO',iwso)
       write(*,*) 'IWU RWU WSO ',iiwu,irwu,iwso
       call SET(fl,fr,fb,ft,ul,ur,ub,ut,ll)
       RETURN
       end

        subroutine lbar2(llb,n1)
          dimension iasf(13)
           integer ll, n1, n2
          real fl,fr,ft,fb,ul,ur,ut,ub
        character*4 llb(16)
             integer lnd(16),lsha(16)
            common /shad/ lnd

          data iasf/13*1/
        call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
        call gsclip(1)
            call gsasf(iasf)
          call gsfais(1)
        call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c       call lblbar(1,0.90,1.,0.1,0.5,n1,0.4,1.,lnd,0,llb,n2,1)
cc       call lblbar(1,0.90,0.96,0.1,0.5,n1,0.6,1.,lnd,0,llb,n2,0)
       rd=float(n1-1)*0.08+0.1
       call lblbar(1,0.90,0.96,rd,rd+0.08,1,0.6,1.,lnd,0,llb,2,0)
c        call lblbar(0,0.1,0.9,.75,0.8,n1,1.,0.4,lnd,0,llb,n2,1)

       call plotif(0.,0.,2)
c       call gsplci(16)
c       call gstxci(16)

       call setusv("LW",3000)
       dy=(0.5-0.1)/float(n1)
       do i=1,n2
       yy=0.1+float(i-1)*dy
c       call plchhq(0.975,yy,llb(i),0.01,0.,0.)
       enddo
        call set(fl,fr,fb,ft,ul,ur,ub,ut,ll)
        return
       end

      subroutine vector2(uu,vv,ix,iy,dummy,iclr,ulen,ivint,ref_sd,rvt,
     +                   ras)
C     usage:  call vector2(uu,vv,ix,iy,dummy,16,0.05,1,10.,3.,0.015)
C     iclr: vector color (follow GSCR setting)
C     ulen:vector reference length ! 0.05  
C     ivint: plot vector interval  ! 1
C     ref_sd:reference speed       ! 10. 
C     rvt: vector thickness        ! 3. 
C     ras: arrow size              ! 0.015
      real uu(ix,iy),vv(ix,iy)
      real px(5),py(5)
      character ca2*2

      call gsplci(iclr)
      IDM=0
      RDM=0.0
      call vvseti('VPO',1)
      call vvseti('SET',0)
      call vvsetr('XC1',1.)
      call vvsetr('XCM',float(ix))
      call vvsetr('YC1',1.)
      call vvsetr('YCN',float(iy))
      call vvseti('SVF',3)
      call vvsetr('USV',dummy)
      call vvsetr('VSV',dummy)

      call vvseti('XIN',ivint)     ! plot vector interval
c      call vvseti('XIN',3)
      call vvseti('YIN',ivint)     ! plot vector interval
      call vvsetc('MNT',' ')
      call vvsetc('MXT',' ')
      call vvinit(uu,ix,vv,ix,RDM,IDM,ix,iy,RDM,IDM)
      call vvsetr('AMN',ras)   ! arrow size
      call vvsetr('LWD',rvt)   ! vector thickness 
c      call vvgetr('DMX',DMX)  ! NDC Maximum Vector Size
      call vvgetr('vmx',vmx)
      call getset(vl,vr,vb,vt,ul,ur,ub,ut,ll)

      dmx=(ulen/ref_sd)*vmx
      vrl=DMX/(vr-vl) 
      call vvsetr('VRL',vrl) ! Vector Reference Length
      call vvectr(uu,vv,RDM,IDM,IDM,RDM)

c      call vvgetr('dmx',dmx)
      call vvseti('vpo',1)

c      ulen=(dmx/vmx)*speed
c      write(*,*)ulen,dmx,vmx

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)

      px(1)=0.88-ulen-0.02
      py(1)=0.11
      px(2)=0.89
      py(2)=0.11
      px(3)=0.89
      py(3)=0.16
      px(4)=0.88-ulen-0.02
      py(4)=0.16
      px(5)=px(1)
      py(5)=py(1)
      call gsfaci(0)
      call gfa(5,px,py)
      call curve(px,py,5)

      xb=0.87-ulen
      xe=0.87
      yb=0.125
      ye=yb

      write(ca2,'(i2)')int(ref_sd)

      lw=2000.*(rvt/1.5)
      call setusv("LW",lw)
      call plchhq( px(1)+0.005 ,yb+0.02 ,ca2//' m s:S:-1',0.008,0.,-1.)
      call vvdraw(xb,yb,xe,ye,ulen,' ',0,0,0,0) ! plot vector reference 

c      xb1=0.17-ulen
c      yb1=0.83
c      xe1=xb1
c      ye1=0.83+ulen
c      call vvdraw(xb1,yb1,xe1,ye1,ulen,' ',0,0,0,0)

      call set(vl,vr,vb,vt,ul,ur,ub,ut,ll)



      return
      end


      subroutine check_ca(i,ii,ist,ied)
c       write(ca4,'(I4)')iix
c       call check_ca(iix,4,ist,ied)
c       call plchhq(xx,0.08,ca4(ist:ied),0.008,0.,0.)
      character*100 ca1
c      character ca2
c     if(ca2.eq.'f' .or. ca2.eq. 'F')
      write(ca1,'(I100)')i 
c      endif

      ied=ii 
      do k=100,100-ii,-1
        if(ca1(k:k).eq.' ')then
        ist=ii+k-100+1
        goto 44
        endif 
      enddo
 44   continue
c      write(*,*)ca1,ist,ied

      return
      end

      subroutine check_ca2(ri,ii,ist,ied)
c       write(ca4,'(I4)')iix
c       call check_ca(iix,4,ist,ied)
c       call plchhq(xx,0.08,ca4(ist:ied),0.008,0.,0.)
      character*100 ca1
c      character ca2
c     if(ca2.eq.'f' .or. ca2.eq. 'F')
      write(ca1,'(f100.1)')ri
c      endif

      ied=ii
      do k=100,100-ii,-1
        if(ca1(k:k).eq.' ')then
        ist=ii+k-100+1
        goto 44
        endif
      enddo
 44   continue
c      write(*,*)ca1,ist,ied

      return
      end

      subroutine p4_interpo(dat,nx,ny,rx,ry,dummy,rout)
      real dat(nx,ny)
      real rx,ry,dummy,rout
      real p1,p2,p3,p4
      ix=rx
      iy=ry

      if(ix.le.0 .or. ix.gt.nx-1 .or.
     +   iy.le.0 .or. iy.gt.ny-1)then
      rout=dummy
      if(ix.eq.nx .and. iy.le.ny)rout=dat(ix,iy)
      if(ix.le.nx .and. iy.eq.ny)rout=dat(ix,iy)

      return
      endif


      p1=dat(ix,iy)
      p2=dat(ix+1,iy)
      p3=dat(ix,iy+1)
      p4=dat(ix+1,iy+1)

      if(p1.eq.dummy.and.p2.ne.dummy)p1=p2
      if(p2.eq.dummy.and.p1.ne.dummy)p2=p1
      v1=p1+(p2-p1)*(rx-float(ix))
      if(p3.eq.dummy.and.p4.ne.dummy)p3=p4
      if(p4.eq.dummy.and.p3.ne.dummy)p4=p3
      v2=p3+(p4-p3)*(rx-float(ix))
      if(v1.eq.dummy.and.v2.ne.dummy)v1=v2
      if(v2.eq.dummy.and.v1.ne.dummy)v2=v1
      rout=v1+(v2-v1)*(ry-float(iy))

      return
      end




      subroutine uv2wdws(u,v,wd,ws)
      real u,v,wd,ws
      pi=4.*atan(1.)
      ws=sqrt(u*u+v*v)
      wd=180.*atan2(-u,-v)/pi
      if(wd.eq.360.)wd=0.
      if(wd.lt.0)wd=wd+360.

      return
      end

      subroutine wdws2uv(wd,ws,u,v)
      real wd,ws,v,u
      pi=4.*atan(1.)
      ctr=pi/180.
      u=ws*cos((270.-wd)*ctr)
      v=ws*sin((270.-wd)*ctr)
      return
      end



      Function Vt(Z,H)

 
C This function computes mean terminal velocity from the reflectivity
C according to Paul Willis' 2-parameter gamma distribution and
C a snow relationship developed by Heymsfield (1978).
C Reflectivity, in terms of dBZ, must by passed to this routine.
C HEIGHT MUST BE IN KM
 
      If (z .gt. 64.0) z = 64.0
      Zz = 10.0**(z/10.0)  ! Reflectivity in mm**6 m**3

C density correction term (rhoo/rho)*0.4 [rho(Z)=rhoo exp-(z/H), where
C  H is the scale height (9.58125) from Gray's inner 2 deg composite]

      Dcor = Exp(0.4*h/9.58125)

C The snow relationship (Atlas et al. 1973) - VT=0.817*Z**0.063  (m/s)

      Vts = -Dcor * (0.817*Zz**0.063)

C The rain relationship -- from Willis' two-parameter analytical
C gamma distribution
C     TERM1=7.331/ZZ**0.010022
C     TERM2=0.14034*ZZ**0.095238
C     VTR=-DCOR * (5.5011E+09/(TERM1+TERM2)**10.5)
C The rain relationship (Joss and Waldvogel,1971)  VT=2.6*Z**.107 (m/s)

      VTR=-DCOR * (2.6*ZZ**.107)

C test if height is in the transition region between SNOW and RAIN
C  defined as 4.7 km < H < 7 km
C  if in the transition region do a linear weight of VTR and VTS

      VT=VTR*(7.-H)/2.3 + VTS*(H-4.7)/2.3
      IF(H.LT.4.7) VT=VTR
      IF(H.GT.7.0) VT=VTS
      Return
      End

      Function variance(var,nxmax,nymax,dummy)
       real var(nxmax,nymax)
       pmean=0.
       ncount=0
       do j=1,nymax
        do i=1,nxmax
         if(var(i,j).ne.dummy)then
          pmean=pmean+var(i,j)
          ncount=ncount+1
         endif
        enddo
       enddo
       if(ncount.ne.0)pmean=pmean/float(ncount)

        pp=0.
       do j=1,nymax
        do i=1,nxmax
         if(var(i,j).ne.dummy)then
          pp=pp+(var(i,j)-pmean)**2
         endif
        enddo
       enddo
        variance=pp/float(ncount)
      Return 
      END

      subroutine axis(nx,nxdi,ccx,ny,nydi,ccy)
c       divide y axis into nx,
c       divide x axis into ny
c      ccx xlabel,ccy ylabel
c      if nxdi,nydi equal 1, labels are below, left
c      if nxdi,nydi equal 2,labels are top, right
      integer n
      character cc*100
      character ccx*(*),ccy*(*)
      real pp,ave_h
      nxdiv=nx
      nydiv=ny
      call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      if(nx.LE.0.and.ny.LE.0)goto 888
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)

      if(nx.GE.1)then
      dx=(fr-fl)/float(nxdiv)
      c_cx=0.001/(fr-fl)
      ave_w=(fr+fl)/2.
c-------------------------------------------------------
      do j=1,(nxdiv+1)
       ddx=ul+(j-1)*((ur-ul)/float(nxdiv))
       write(cc,'(i10)')int(ddx)
       rx=fl+float(j-1)*dx
       if(nxdi.EQ.1)then
       call plchhq(rx,(fb-0.02),trim(adjustl(cc)),0.01,0.,0.)
       if(j.eq.1)call plchhq(ave_w,(fb-0.07),ccx,0.015125-c_cx,0.,0.)
       elseif(nxdi.EQ.2)then
       call plchhq(rx,(ft+0.02),trim(adjustl(cc)),0.01,0.,0.)
       if(j.eq.1)call plchhq(ave_w,(ft+0.07),ccx,0.015125-c_cx,0.,0.)
       endif
      enddo
      endif
c--------------------------------------------------------
      if(ny.GE.1)then
      dy=(ft-fb)/float(nydiv)
      c_cy=0.001/(ft-fb)
      ave_h=(ft+fb)/2.
      do j=1,(nydiv+1)
        ddy=ub+(j-1)*((ut-ub)/float(nydiv))
        write(cc,'(f10.1)')ddy
        ry=fb+float(j-1)*dy
        if(nydi.EQ.1)then
        call plchhq((fl-0.005),ry,trim(adjustl(cc)),0.01,0.,1.)
        if(j.eq.1)call plchhq((fl-0.08),ave_h,ccy,0.015125-c_cy,90.,0.)
        elseif(nydi.EQ.2)then
        call plchhq((fr+0.01),ry,trim(adjustl(cc)),0.013,0.,-1.)
        if(j.eq.1)call plchhq((fr+0.08),ave_h,ccy,0.015125-c_cy,90.,0.)
        endif
        enddo
      endif
c--------------------------------------------------------
        call set(fl,fr,fb,ft,ul,ur,ub,ut,ll)
        call perim(abs(nxdiv),2,abs(nydiv),2)
888     continue
        return
        end

