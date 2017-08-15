       Program plot_dual_doppler_outdata
       integer ixmax,iymax
       character filein*100       
C------------------------------------------------------------

       write(*,*)'input CEDREAD outdata ='
       read(*,'(A100)')filein
       open(10,file=filein,status='old')
       read(10,*)ixmax,iymax,level,ifield
       close(10)
       call plot_dual(ixmax,iymax,filein)
       end
       
       subroutine plot_dual(ixmax,iymax,filein)   
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
       integer domain
       integer iasf(13),ity(21)
       character llb(16)*2,ca7*14,ca8*14,ca9*14,date*12
       character ca2*2,ca4*4,ca1*1,ca3*3,ca5*5,case_name*8
       character filein*100,cyy*4,cmo*2,cdd*2,chh*2,cmm*2,clvl*5
       common /pattern/ity
       data iasf/13*1/
       data ity/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21/
       real external Vt  ! terminal velocity calc
       real external variance
       real*8 xdis_rd(2),ydis_rd(2),zdis_rd(2),R(2),vr_vri(2)
       real*8 Sxx,Syy,Sxy,Sxz,Syz,tmp,vtmp,wtmp

C-------------- pattern center -----------------------
       ca7='$$$$$$$$$$$$$$'
       ca8="$'$'$'$'$'$'$'"
       ca9="$$$$'$'$$$$'$'"
       pi=4.*atan(1.)
       ctr=pi/180.
       dummy=-1000.

       open(40,file='ty_info4lwcheng2_new.txt',status='unknown')
       read(40,*)
       do i=1,100
       read(40,*,end=999)case_name,date,ty_theta,ty_speed,domain
       if(case_name(1:8).EQ.filein(1:8))then
       go to 999
       endif
       enddo
999    continue
       print *,'The rainband speed,dir:',ty_speed,ty_theta
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
        west=20.
        east=100.
        south=-120.
        rnorth=-10.
       elseif(domain.EQ.3)then
       ox=299701.532  ! GI 
       oy=2505893.566
        west=25.
        east=85.
        south=30.
        rnorth=100.
       elseif(domain.EQ.4)then
       ox=234181.303  ! kt 
       oy=2422789.221
        west=60.
        east=130.
        south=-20.
        rnorth=70.
       elseif(domain.EQ.5)then
       ox=156399.173  ! cg
       oy=2560832.205
        west=-50.
        east=10.
        south=-130.
        rnorth=-70.
       elseif(domain.EQ.6)then
       ox=156399.173  ! cg
       oy=2560832.205
        west=-110.
        east=-30.
        south=-60.
        rnorth=20.
       elseif(domain.EQ.7)then
       ox=156399.173  ! cg
       oy=2560832.205
        west=-70.
        east=0.
        south=60.
        rnorth=130.
       elseif(domain.EQ.8)then
       ox=270176.646 ! cks  
       oy=2773555.278
        west=-110.
        east=-40.
        south=-30.
        rnorth=30.
       endif
C------------------------------
       cyy(1:4)=date(1:4)
       cmo(1:2)=date(5:6)
       cdd(1:2)=date(7:8)
       chh(1:2)=date(9:10)
       cmm(1:2)=date(11:12)
c-----------------------------------------
       write(*,*)'input plot mode =[1:vertical wind 2:dBZ&wind]'
       read(*,*)pmode

       write(*,*)'input LEVEL(m)='
       read(*,*)nlevel

       write(*,*)'you input level =',nlevel,' m'
       nlevel=nlevel/500
       open(30,file=filein,status='old')

       do k=1,nlevel
       read(30,*)nxmax,nymax,level,ifield
         if(nxmax.ne.ixmax  .or. nymax.ne.iymax)then
         write(*,*)'ixmax or iymax WRONG !!!'
         write(*,*)'modify ixmax and iymax to',nxmax,nymax
         call exit(1)
         endif
C----------------------- 3 
       do j=1,nymax
       do i=1,nxmax
       read(30,*)dbz(i,j)
       enddo
       enddo

C----------------------- 8
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)
       enddo
       enddo
C----------------------- 9
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)vr_radar1(i,j)
       enddo
       enddo
C----------------------- 10
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)
       enddo
       enddo
C----------------------- 11 
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)vr_radar2(i,j)
       enddo
       enddo

C----------------------- 12
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)wfv(i,j)
       enddo
       enddo

C----------------------- 14
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)
       enddo
       enddo

C----------------------- 16
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)uf(i,j)
       enddo
       enddo

C----------------------- 17
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)vf(i,j)
       enddo
       enddo


c--------subtract the speed of typhoon rainband-------
        ty_speed=0.0 !testing
        ty_u=ty_speed*cos(ty_theta*pi/180.)
        ty_v=ty_speed*sin(ty_theta*pi/180.)
       print *,ty_speed,ty_u,ty_v
       do j=1,nymax
       do i=1,nxmax
       if(uf(i,j).ne.dummy .and. vf(i,j).ne. dummy)then 
       uf(i,j)=uf(i,j)-ty_u
       vf(i,j)=vf(i,j)-ty_v
       endif
       enddo
       enddo


       enddo
       close(30)
c=====calculate vt using height and dBZ
      pheight=float(nlevel*500) !the point height
      vtv=dummy
      do j=1,nymax
       do i=1,nxmax
        if(dbz(i,j).ne.dummy)then
        vtv(i,j)=Vt(dbz(i,j),(pheight/1000.))
        endif
       enddo
      enddo
c     calculate variance of vr,w,vt
      vr_vri(1)=variance(vr_radar1,nxmax,nymax,dummy)
      vr_vri(2)=variance(vr_radar2,nxmax,nymax,dummy)
      vtv_vri=variance(vtv,nxmax,nymax,dummy)
      w_vri=variance(wfv,nxmax,nymax,dummy)
c=====calculate variance==============
      
c     the point where the radial wind is observed
      ptx=-30.
      pty=50.
      pox=ox + (ptx*1000.)
      poy=oy + (pty*1000.)
c     calc the distance x y z
c     xdis_rd == x-xi , ydis_rd == y-yi , zdir_rd == z-zi
      xdis_rd(1) = (pox-ox)  ! distance btw rd1 and pt
      ydis_rd(1) = (poy-oy) 
      zdis_rd(1) = (pheight-oz)
      xdis_rd(2) = (pox-ox2) ! distance btw rd2 and pt
      ydis_rd(2) = (poy-oy2)
      zdis_rd(2) = (pheight-oz2)
c      print *,xdis_rd,ydis_rd,zdis_rd,'xyz'
c     R == Ri  R(1)=radar1, R(2)=radar2
      do i=1,2
       R(i) = (xdis_rd(i)**2+ydis_rd(i)**2+zdis_rd(i)**2)**0.5
      enddo
      print *,R,'~~~~'
c     calc Sxx, Syy, Sxy, Sxz, Syz
      Sxx = 0.
      Syy = 0.
      Sxy = 0.
      Sxz = 0.
      Syz = 0.
      do i=1,2
       Sxx = Sxx + xdis_rd(i)**2
       Syy = Syy + ydis_rd(i)**2
       Sxy = Sxy + xdis_rd(i)*ydis_rd(i)
       Sxz = Sxz + xdis_rd(i)*zdis_rd(i)
       Syz = Syz + ydis_rd(i)*zdis_rd(i)
      enddo
       print *,Sxx,Syy,Sxy,Sxz,Syz
c=====calculate variance of u,v,w
       utmp=0.
       vtmp=0.
       wtmp=0.
      do i=1,2
       utmp = (utmp +
     &   vr_vri(i)*(R(i)*xdis_rd(i)*Syy-R(i)*ydis_rd(i)*Sxy)**2)
     &   /((Sxx*Syy-(Sxy)**2)**2)
       vtmp = (vtmp +
     &   vr_vri(i)*(R(i)*ydis_rd(i)*Sxx-R(i)*xdis_rd(i)*Sxy)**2)
     &  /((Sxx*Syy-(Sxy)**2)**2)
      print *,utmp,vtmp,'utmp vtmp'
      enddo
      u_vari = utmp+((w_vri+vtv_vri)*((Sxy*Syz-Syy*Sxz))**2
     &        /((Sxx*Syy-(Sxy)**2))**2)
      v_vari = vtmp+((w_vri+vtv_vri)*((Sxy*Sxz-Sxx*Syz))**2
     &        /((Sxx*Syy-(Sxy)**2))**2)
      print *,u_vari,v_vari,'u_vari v_vari'
c=====end of codes====================




       call opngks
       call defclr16
       call gsplci(16)
       call setusv('LW',2000)
c       call gsfais(1)
c       call gstxci(16)
       call gsclip(1)

        do i=1,16
        wfv1_clv(i)=0.+0.5*float(i-1)
        wfv2_clv(i)=0.-0.5*float(i)
        wfv_clv(i)=-3.5+0.5*float(i-1)
        enddo
        wfv_clv(1)=-30.
        wfv_clv(15)=30.
        wfv_clv(16)=50.

        do i=1,21
        dbz_clv(i)=-5.+5.*float(i-1)
        write(llb(i),'(i2)')int(dbz_clv(i))
        enddo
c        dbz_clv(13)=70
        dbz_clv(1)=-50.
        llb(1)='  '
        llb(15)='  '

       
       if(iymax.gt.ixmax)then
       dd=0.8*(float(ixmax-1)/float(iymax-1))

       call set((0.5-dd/2.),(0.5+dd/2.),0.1,0.9,
     +           1.,float(nxmax),1.,float(nymax),1)

       elseif(iymax.eq.ixmax)then
       dd=0.8*(float(ixmax-1)/float(iymax-1))
       call set(0.1,0.1+dd,0.1,0.9,
     +           1.,float(nxmax),1.,float(nymax),1) 
       else 
       dd=0.8*(float(iymax-1)/float(ixmax-1))
       call set(0.1,0.9,0.1,0.1+dd,
     +           1.,float(nxmax),1.,float(nymax),1)
       endif

        call lbar(llb,12,13,pmode)

       if(pmode .eq.1 )then

       call color(dbz,ixmax,iymax,dbz_clv,dummy)
       call contour(wfv,ixmax,iymax,dummy,-6.,1,1.5,ca8,16,4)
       call contour(wfv,ixmax,iymax,dummy,-5.,1,1.5,ca8,16,4)
       call contour(wfv,ixmax,iymax,dummy,-4.,1,1.5,ca8,16,4)
       call contour(wfv,ixmax,iymax,dummy,-3.,1,1.5,ca8,16,4)
       call contour(wfv,ixmax,iymax,dummy,-2.,1,1.5,ca8,16,4)
       call contour(wfv,ixmax,iymax,dummy,-1.,1,1.5,ca8,16,4)
       call contour(wfv,ixmax,iymax,dummy,0.,1,1.5,ca9,16,4)
       call contour(wfv,ixmax,iymax,dummy,1.,1,1.5,ca7,16,4)
       call contour(wfv,ixmax,iymax,dummy,2.,1,1.5,ca7,16,4)
       call contour(wfv,ixmax,iymax,dummy,3.,1,1.5,ca7,16,4)
       call contour(wfv,ixmax,iymax,dummy,4.,1,1.5,ca7,16,4)
       call contour(wfv,ixmax,iymax,dummy,5.,1,1.5,ca7,16,4)
       call contour(wfv,ixmax,iymax,dummy,6.,1,1.5,ca7,16,4)

       elseif(pmode .eq. 2)then

       call color(dbz,ixmax,iymax,dbz_clv,dummy)
       call vector2(uf,vf,ixmax,iymax,dummy,16,0.05,3,30.,3.,0.015)
c       call vector2(uf,vf,ixmax,iymax,dummy,16,0.07,3,30.,3.,0.015)
       ref_w=maxval(uf)
       ref_w2=maxval(vf)
       ref_w3=max(ref_w,ref_w2)
       ref_w3=ref_w3
       print *,ref_w3
c       call vector2(uf,vf,ixmax,iymax,dummy,16,0.05,3,ref_w3,3.,0.015) 
       endif



C--------------------- to plot Taiwan outline and topo-----------------
      call shader_topo(ox,oy,east,west,south,rnorth,16)
C--------------------- to plot Taiwan outline and topo done-----------------

      ix=nint(east-west)/10
      iy=nint(rnorth-south)/10

       call gsplci(16)
       call setusv('LW',2000)
       call perim(ix,5,iy,5)

       call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
       call set(0.,1.,0.,1.,0.,1.,0.,1.,1)

 
       if(iymax.ge.ixmax)then
       dd1=dd/float(ix)
       else
       dd1=0.8/float(ix)
       endif
       do i=1,ix+1
c       ixx=west+float(i-1)*10.
       ixx=0.+float(i-1)*10.
       write(ca4,'(I4)')ixx
       if(iymax.le.ixmax)then
       xx=0.1+float(i-1)*dd1
       else
       xx=(0.5-dd/2.)+float(i-1)*dd1
       endif
       call check_ca(ixx,4,ist,ied)
       call plchhq(xx,0.08,ca4(ist:ied),0.015,0.,0.)
       enddo

       if(iymax.ge.ixmax)then
       dd1=0.8/float(iy)
       else
       dd1=dd/float(iy)
       endif

       do i=1,iy+1
c       iyy=south+float(i-1)*10.
       iyy=0.+float(i-1)*10.
       write(ca4,'(I4)')iyy
       yy=0.1+float(i-1)*dd1
       if(iymax.le.ixmax)then
       call plchhq(0.06,yy,ca4,0.015,0.,0.)
       else
       call plchhq((0.46-dd/2.),yy,ca4,0.015,0.,0.)
       endif
       enddo
       

        
       write(ca5,'(f5.2)')float(nlevel)*0.5
       call plchhq(0.72,ft+0.015,ca5,0.015,0.,0.)
       call plchhq(0.84,ft+0.015,'km CAPPI',0.015,0.,0.)
       call plchhq((fl+fr)/2.,0.03,'X (km)',0.018,0.,0.)
       if(iymax.le.ixmax)then
       call plchhq(0.025,(ft+fb)/2.,'Y (km)',0.018,90.,0.)
       else
       call plchhq((0.425-dd/2.),(ft+fb)/2.,'Y (km)',0.018,90.,0.)
       endif
C------------------- plot time -------------------------
       if(1.eq.1)then
       call plchhq(0.62,0.03,cyy,0.015,0.,0.)
       call plchhq(0.700,0.03,cmo,0.015,0.,0.)
       call plchhq(0.725,0.03,'/',0.015,0.,0.)
       call plchhq(0.750,0.03,cdd,0.015,0.,0.)
       call plchhq(0.795,0.03,chh,0.015,0.,0.)
       call plchhq(0.815,0.03,':PRL:0',0.015,0.,0.)
       call plchhq(0.835,0.03,cmm,0.015,0.,0.)
       call plchhq(0.884,0.03,'UTC',0.015,0.,0.)
       endif



C------------------- RHI start
       write(*,*)'RHI mode ? [yes:1 ; no:2]'
       read(*,*)pmode
       if(pmode.eq.1)then
       write(*,*)'input rx and ry (km)'
       read(*,*)px,py
       write(*,*)'input length(km) and height(km)'
       read(*,*)leng,km
       write(*,*)'decide theta automatically?[yes:1; no:other number]'
       read(*,*)pmode2
C----------find the nearest boundary to cross------------------
       if(pmode2.eq.1)then
c----------------------------------------
       else
       write(*,*)'Please enter the theta'
       read(*,*)theta
       theta=theta+180.
       endif
c-----------set----------------------       
       nnz=20
       hor_space=1.
       ver_space=0.5
C--------------------------
       ixt=(ixmax-1)*hor_space*2.5+1
       iyt=(iymax-1)*hor_space*2.5+1

       leng1=leng/hor_space+1
       leng2=leng*2.5+1
       km1=km/ver_space+1

      call set(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      call setusv("LW",9000)
      call setusv("LW",2000)
      call plot_line(px,py,hor_space,theta,leng1,px4,py4)
      call frame
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
       write(*,*)'only plot the second figure ?'
     +         //'[yes:1 no:any other number]'
       read(*,*)ipmode
      if(ipmode.eq.1)then
      call clsgks
      call opngks
      call defclr16
      call setusv("LW",2000)
      call gsplci(16)
      call gstxci(16)
      endif
C-----------------------------------------

      leng1=(leng1-1)*2+1

c      call rhi(filein,ox,oy,ixmax,iymax,nnz,hor_space,ver_space,
c     +           rnorth,south,east,west,ixt,iyt,
c     +           nlevel,px,py,theta,leng1,leng2,km1,px4,py4,
c     +           ty_speed,ty_theta,date)
       endif
C------------------- RHI end

       call frame
       call clsgks

       stop
       end

       subroutine plot_line(px,py,hor_space,dir,leng1,px4,py4)
       pi=4.*atan(1.)
       ctr=pi/180.
c       px1=px/hor_space+1.
c       py1=py/hor_space+1.
        px1=px/hor_space
        py1=py/hor_space

       do ii=1,1
       px2=px1+float(ii-1)*sin((dir-90.)*ctr)
       py2=py1+float(ii-1)*cos((dir-90.)*ctr)
       
       do jj=1,leng1,leng1-1
       px3=px2+float(jj-1)*sin(dir*ctr)
       py3=py2+float(jj-1)*cos(dir*ctr)
       
       px4=px2-float(jj-1)*sin(dir*ctr)
       py4=py2-float(jj-1)*cos(dir*ctr)
       enddo
       enddo
       call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
       dd2=5.
       ffl=(fr-fl)*dd2/(ur-ul)
       ffb=(ft-fb)*dd2/(ut-ub)
       call set((fl-ffl),(fr+ffl),(fb-ffb),(ft+ffb),
     &           (ul-dd2),(ur+dd2),(ub-dd2),(ut+dd2),ll)
       call setusv("LW",9000)
       call line(px4,py4,px3,py3)
       call gsplci(16)
       call setusv("LW",9000)
       call plchhq((px4-3.*sin(dir*ctr)),(py4-3.*cos(dir*ctr))
     &              ,'A',0.02,0.,0.)
       call plchhq((px3+3.*sin(dir*ctr)),(py3+3.*cos(dir*ctr))
     &               ,"A'",0.02,0.,0.)
       call gsplci(3)
       call setusv("LW",20000)
       call point(px1,py1)
       call gsplci(16)
c       print *,((px4-px3)**2+(py4-py3)**2)**(0.5)
c       print *,((px2-px3)**2+(py2-py3)**2)**(0.5)
c       call gsplci(2)
c       call line(px2,py2,px3,py3)
       call setusv("LW",2000)
       call set(fl,fr,fb,ft,ul,ur,ub,ut,ll)

       return
       end

       subroutine rhi(filein,ox,oy,nnx,nny,nnz,hor_space,ver_space,
     +                 north,south,east,west,ix,iy,
     +                 nlevel,px,py,dir,leng1,leng2,km1,px4,py4,
     +                 ty_speed,ty_theta,date)
       real grid3dz(nnx,nny,nnz)
       real grid3du(nnx,nny,nnz)
       real grid3dv(nnx,nny,nnz)
       real grid3dw(nnx,nny,nnz)
       real topo(ix,iy),clv_height(16)
       real clv_ref(16)
       real clv_ws(16)
       character filein*100
       parameter(dummy=-1000.)
       real csz(leng1,nnz+1),topo2(leng2+2),topox(leng2+2)
       real csvp(leng1,nnz+1),csvc(leng1,nnz+1)
       real csw(leng1,nnz+1)
c       real ws(leng1,nnz)
       real clv_pp(21)
       real x(10000),rx(10000),ry(10000)
       real*8 gaea_x,gaea_y
       real ty_speed,ty_u,ty_v
       integer distance,iyd,ih,im,topo_dis
       integer iasf(13),ity(21),nlevel
       character llb(21)*2,ca7*14,ca8*14,ca9*14,ca_date*6,date*12
       character llb2(16)*4
       character cyy*4,cmo*2,cdd*2,chh*2,cmm*2,clvl*5
       character ca4*4
       character filename*100,ca_2*2,ca2*2
       data iasf/13*1/
       data ity/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21/
       common /pattern/ity
       integer lsha(16)

        pi=4.*atan(1.)
        ctr=pi/180.

       open(30,file=filein,status='old')

        
       do k=1,nnz
       read(30,*)nxmax,nymax,level,ifield
         if(nxmax.ne.nnx  .or. nymax.ne.nny)then
         write(*,*)'ixmax or iymax WRONG !!!'
         write(*,*)'modify ixmax and iymax to',nxmax,nymax
         call exit(1)
         endif
C----------------------- 3 
       do j=1,nymax
       do i=1,nxmax
       read(30,*)grid3dz(i,j,k)
       enddo
       enddo

C----------------------- 8
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)
       enddo
       enddo
C----------------------- 9
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)
       enddo
       enddo
C----------------------- 10
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)
       enddo
       enddo
C----------------------- 11 
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)
       enddo
       enddo

C----------------------- 12
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)grid3dw(i,j,k)
       enddo
       enddo

C----------------------- 14
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)
       enddo
       enddo

C----------------------- 16
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)grid3du(i,j,k)
       enddo
       enddo

C----------------------- 17
       read(30,*)nxmax,nymax,level,ifield
       do j=1,nymax
       do i=1,nxmax
       read(30,*)grid3dv(i,j,k)
       enddo
       enddo

c--------subtract the speed of typhoon rainband for cross section data-------
       ty_u=ty_speed*cos(ty_theta*pi/180.)
       ty_v=ty_speed*sin(ty_theta*pi/180.)
       print *,ty_speed,ty_theta
       do j=1,nymax 
       do i=1,nxmax 
       if(grid3du(i,j,k).ne.dummy .and. grid3dv(i,j,k).ne. dummy)then
       grid3du(i,j,k)=grid3du(i,j,k)-ty_u
       grid3dv(i,j,k)=grid3dv(i,j,k)-ty_v 
       endif
       enddo
       enddo


       enddo


C------------------ read 3D data ok

C------------------ cross section data -------------
       px1=px/hor_space+1.
       py1=py/hor_space+1.

       csz=dummy
       csvp=dummy
       csvc=dummy
       csw=dummy
       
       do k=1,nnz
       do j=1,leng1
      px2=px4+float(j-1)*sin((dir)*ctr)
      py2=py4+float(j-1)*cos((dir)*ctr)

c       write(*,*)px2,py2

      call p4_interpo(grid3dz(1:nnx,1:nny,k:k),nnx,nny,
     + px2,py2,dummy,csz(j,k+1))
      call p4_interpo(grid3dw(1:nnx,1:nny,k:k),nnx,nny,
     + px2,py2,dummy,csw(j,k+1))

      call p4_interpo(grid3du(1:nnx,1:nny,k:k),nnx,nny,
     + px2,py2,dummy,uu)
      call p4_interpo(grid3dv(1:nnx,1:nny,k:k),nnx,nny,
     + px2,py2,dummy,vv)
       if(uu.ne.dummy.and.vv.ne.dummy)then
       call uv2wdws(uu,vv,wd,ws)
       diff=wd-dir
       csvp(j,k+1)=-ws*cos(diff*ctr)
       csvc(j,k+1)= ws*sin(diff*ctr)
       endif

       enddo
       enddo
C------------------ cross section data done-------------
C------------------ cross section topo data -------------
       call get_topo_400(ox,oy,west,south,topo,ix,iy,2)
       px1=px*2.5+1.
       py1=py*2.5+1.
       do j=1,leng2
       px2=px1+float(j-1)*sin((dir)*ctr)
       py2=py1+float(j-1)*cos((dir)*ctr)
       call p4_interpo(topo,ix,iy,px2,py2,dummy,topo2(j))
       enddo
C------------------ cross section topo data done-------------
C--------------- plot
      call setusv("LW",2000)

      do i=1,21
      clv_pp(i)=-5.+float(i-1)*5.
      write(llb(i),'(i2)')int(clv_pp(i))
      enddo
      clv_pp(1)=-100.
      llb(1)='   '
       ref_wind=0.
       do i=1,leng1
       do j=1,km1
       if(csw(i,j).ne.dummy)then
c       csw(i,j)=csw(i,j)*2.
       endif
       if(csvp(i,j).ne.dummy)then
        if(abs(csvp(i,j)).GT.ref_wind)then
        ref_wind=abs(csvp(i,j))
        endif
c       csvp(i,j)=csvp(i,j)*2.
       endif
       enddo
       enddo
       call set(0.1,0.9,0.1,0.9,1.,float(leng1),1.,float(km1),1)
       call color(csz(1:leng1,1:km1),leng1,km1,clv_pp,dummy)
       call lbar(llb,12,13,2)
       call vector2(csvp(1:leng1,1:km1),csw(1:leng1,1:km1),leng1,km1,
     +              dummy,16,0.05,1,(ref_wind/1.5),3.,0.015)
c       call vector2(csvp(1:leng1,1:km1),csw(1:leng1,1:km1),leng1,km1,
c     +              dummy,16,0.05,2,4.,3.,0.015)
       print *,'the reference wind is:',(ref_wind/1.5)
c---------------------------------------------------------
       
       open(15,file=filein(1:8)//'_dBZ.txt')
       open(16,file=filein(1:8)//'_hori.txt')
       open(17,file=filein(1:8)//'_verti.txt')
       do ii=1,leng1
        do jj=1,km1
       write(15,*)csz(ii,jj)
       write(16,*)csvp(ii,jj)
       write(17,*)csw(ii,jj)
         enddo
       enddo
       close(15)
       close(16)
       close(17)
C--------------- plot topo
       call set(0.1,0.9,0.1,0.9,1.,float(leng2)
     +           ,0.,float(km1-1)*ver_space*1000.,1)

       topox(1)=1.
       do i=2,leng2+1
       topox(i)=float(i-1)
       enddo
       topox(leng2+2)=float(leng2)

       do i=leng2,1,-1
       topo2(i+1)=topo2(i)
       enddo
       topo2(1)=0.
       topo2(leng2+2)=0.

       call GSFAIS(1) !area filling's type [1 or 3]
       call GSFACI(22)  !change area filling's color
c       call GFA(leng2+2,topox,topo2)



C-------------------------------------
      call setusv("LW",2000)
       call perim(1,1,1,1)
C---------------------------  label ------------------------------
       ileng=(leng1-1)*hor_space

       if(mod(ileng,10).eq.0)then
       idx1=ileng/10
       idx2=10
       endif
c       print *,leng1,ileng,idx1


       call perim(idx1,idx2,km1-1,1)
       call set(0.,1.,0.,1.,0.,1.,0.,1.,1)


       dd1=0.8/(idx1)
       do i=1,idx1+1
        iix=0+(i-1)*idx2
        write(ca4,'(i4)')iix
       call check_ca(iix,4,ist,ied)
        xx=0.1+(i-1)*dd1
       call plchhq(xx,0.08,ca4(ist:ied),0.012,0.,0.)
       enddo



       dd1=0.8/float(km1-1)
       do i=1,km1
        rry=0.+float(i-1)*ver_space
        iiy=rry
        if(float(iiy).eq.rry)then
        write(ca4,'(I4)')iiy
        xx=0.1+(i-1)*dd1
        call plchhq(0.066,xx,ca4,0.015,0.,0.)
        endif
       enddo
       call plchhq(0.1,0.04,'Front(A)',0.015,0.,0.)
       call plchhq(0.9,0.04,"Rear(A')",0.015,0.,0.)
C---------------------------  plot time ----------------------------
       cyy(1:4)=date(1:4)
       cmo(1:2)=date(5:6)
       cdd(1:2)=date(7:8)
       chh(1:2)=date(9:10)
       cmm(1:2)=date(11:12)
       if(1.eq.1)then
       call plchhq(0.61,0.03,cmo,0.015,0.,0.)
       call plchhq(0.635,0.03,'/',0.015,0.,0.)
       call plchhq(0.660,0.03,cdd,0.015,0.,0.)
       call plchhq(0.705,0.03,chh,0.015,0.,0.)
       call plchhq(0.725,0.03,':PRL:0',0.015,0.,0.)
       call plchhq(0.745,0.03,cmm,0.015,0.,0.)
       call plchhq(0.794,0.03,'UTC',0.015,0.,0.)
       endif



       call plchhq(0.5,0.03,'X (km)',0.018,0.,0.)
       call plchhq(0.04,0.5,'Height (km)',0.018,90.,0.)

        

C-----------------------------------

      return 
      end



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
       dimension rwrk(40000),iwrk(30000) !by yuku
c       dimension rwrk(160000),iwrk(160000) !by lwcheng
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
c          call CPSETI('CLU - CONTOUR LEVEL USE',3)
         call CPSETI('CLU - CONTOUR LEVEL USE',1)
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
