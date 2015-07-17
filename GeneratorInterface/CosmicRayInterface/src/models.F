      subroutine NumberModel(cmodel,model)
      character cmodel*21
      n=index(cmodel,' ')-1
      if(cmodel(1:n).eq.'epos' )     model=1
      if(cmodel(1:n).eq.'qgsjet' )then
#ifndef __QGSJET01__
         stop'please compile with requested model'
#endif
         model=2
      endif
      if(cmodel(1:n).eq.'gheisha' )then
#ifndef __GHEISHA__
         stop'please compile with requested model'
#endif
         model=3
      endif
      if(cmodel(1:n).eq.'pythia' )then
#ifndef __PYTHIA__
         stop'please compile with requested model'
#endif
         model=4
      endif
      if(cmodel(1:n).eq.'hijing' )then
#ifndef __HIJING__
         stop'please compile with requested model'
#endif
         model=5
      endif
      if(cmodel(1:n).eq.'sibyll' )then
#ifndef __SIBYLL__
         stop'please compile with requested model'
#endif
         model=6
      endif
      if(cmodel(1:n).eq.'IIqgsjet' )then
#ifndef __QGSJETII__
         stop'please compile with requested model'
#endif
         model=7
      endif
      if(cmodel(1:n).eq.'phojet' )then
#ifndef __Phojet__
         stop'please compile with requested model'
#endif
         model=8
      endif
      if(cmodel(1:n).eq.'fluka' )then
         stop'wrong choice!!!!' !    model=9
      endif
      end

      subroutine IniModel(model)
      if(model.eq.2) then
#ifndef __QGSJET01__
         stop'please compile with requested model'
#endif
         call IniQGSjet
      endif
      if(model.eq.3) then
#ifndef __GHEISHA__
         stop'please compile with requested model'
#endif
         call IniGheisha
      endif
      if(model.eq.4) then
#ifndef __PYTHIA__
         stop'please compile with requested model'
#endif
         call IniPythia
      endif
      if(model.eq.5) then
#ifndef __HIJING__
         stop'please compile with requested model'
#endif
         call IniHijing
      endif
      if(model.eq.6) then
#ifndef __SIBYLL__
         stop'please compile with requested model'
#endif
         call IniSibyll
      endif
      if(model.eq.7) then
#ifndef __QGSJETII__
         stop'please compile with requested model'
#endif
         call IniQGSJetII
      endif
      if(model.eq.8) then
#ifndef __Phojet__
         stop'please compile with requested model'
#endif
         call IniPhojet
      endif
c     if(model.eq.9)call IniFluka
      end

      subroutine IniEvtModel
      include 'epos.inc'
      iclegy=1                  !to avoid crash in plots
      if(model.eq.2)then
#ifndef __QGSJET01__
         stop'please compile with requested model'
#endif
         call IniEvtQGS
      endif
      if(model.eq.3)then
#ifndef __GHEISHA__
         stop'please compile with requested model'
#endif
         call IniEvtGhe
      endif
      if(model.eq.4)then
#ifndef __PYTHIA__
         stop'please compile with requested model'
#endif
         engysave=engy
         if(engy.lt.egymin)engy=egymin
         call IniEvtPyt
         engy=engysave
      endif
      if(model.eq.5)then
#ifndef __HIJING__
         stop'please compile with requested model'
#endif
         engysave=engy
         if(engy.lt.egymin)engy=egymin
         call IniEvtHij
         engy=engysave
      endif
      if(model.eq.6)then
#ifndef __SIBYLL__
         stop'please compile with requested model'
#endif
         call IniEvtSib
      endif
      if(model.eq.7)then
#ifndef __QGSJETII__
         stop'please compile with requested model'
#endif
         call IniEvtQGSII
      endif
      if(model.eq.8)then
#ifndef __Phojet__
         stop'please compile with requested model'
#endif
         call IniEvtPho
      endif
c     if(model.eq.9)call IniEvtFlu
      end

      subroutine emsaaaModel(model,id,iret)
      if(model.eq.2) then
#ifndef __QGSJET01__
         stop'please compile with requested model'
#endif
         if(id.eq.0)call IniEvtQGS
         call emsqgs(iret)
      endif
      if(model.eq.3) then
#ifndef __GHEISHA__
         stop'please compile with requested model'
#endif
         call emsghe(iret)
      endif
      if(model.eq.4) then
#ifndef __PYTHIA__
         stop'please compile with requested model'
#endif
         call emspyt(iret,0)
      endif
      if(model.eq.5) then
#ifndef __HIJING__
         stop'please compile with requested model'
#endif
         call emshij(iret)
      endif
      if(model.eq.6) then
#ifndef __SIBYLL__
         stop'please compile with requested model'
#endif
         call emssib(iret)
      endif
      if(model.eq.7) then
#ifndef __QGSJETII__
         stop'please compile with requested model'
#endif
         if(id.eq.0)call IniEvtQGSII
         call emsqgsII(iret)
      endif
      if(model.eq.8) then
#ifndef __Phojet__
         stop'please compile with requested model'
#endif
         call emspho(iret)
      endif
c     if(model.eq.9)call emsflu(iret)
      end


      subroutine crseaaModel(sigt,sigi,sigc,sige)
      include 'epos.inc'
      sigt=0.
      sigc=0.
      sigi=0.
      sige=0.
      if(model.eq.3.or.model.eq.4.or.model.eq.8)return !no AA with Gheisha, Pythia, Phojet
      if(idtarg.eq.0.and.model.ne.9)then
         kmax=3
         if(model.eq.6)kmax=2   !no Argon with SIBYLL
         do k=1,kmax
            matarg=int(airanxs(k))
            call crseaaModel0(xsigt,xsigi,xsigc,xsige)
            sigt=sigt+airwnxs(k)*xsigt
            sigi=sigi+airwnxs(k)*xsigi
            sigc=sigc+airwnxs(k)*xsigc
            sige=sige+airwnxs(k)*xsige
         enddo
      else
         call crseaaModel0(sigt,sigi,sigc,sige)
      endif
      end

      subroutine crseaaModel0(sigt,sigi,sigc,sige)
      include 'epos.inc'

      double precision GTOT,GPROD,GABS,GDD,GQEL,GCOH
      dimension dumdif(3)
      if(model.eq.2)then
#ifndef __QGSJET01__
         stop'please compile with requested model'
#endif
         NITER=5000
         if(idtarg.eq.0)then
            e0=dble(elab)
            icp=idtrafo('nxs','qgs',idproj)
            call xxaini(e0,icp,maproj,matarg)
         endif
         CALL CROSSC(NITER,GTOT,GPROD,GABS,GDD,GQEL,GCOH)
         sigt=sngl(GTOT)
         sigi=sngl(GPROD)
         sigc=sngl(GABS)
         sige=sigt-sigi
      elseif(model.eq.3)then
#ifndef __GHEISHA__
         stop'please compile with requested model'
#endif
         idtar=idtarg
         if(idtarg.eq.0)idtar=1120
         call ghecrse(ekin,idproj,idtar,latarg,matarg,sigi,sige)
         sigt=sigi+sige
         sigc=sigi
      elseif(model.eq.5)then
#ifndef __HIJING__
         stop'please compile with requested model'
#endif
         call hjcrossc(sigi,sigt)
      elseif(model.eq.6.and.maproj.eq.1)then
#ifndef __SIBYLL__
         stop'please compile with requested model'
#endif
         K=1
         if(iclpro.eq.1)then
            K=2
         elseif(iclpro.eq.3)then
            K=3
         endif
         CALL SIB_SIGMA_HP(K,engy,SSIG,dum0,dum1,dumdif,SLOPE,RHO)
         CALL GLAUBER(matarg,SSIG,SLOPE,RHO,sigt,sige,sigqel)
         sigi=sigt-sigqel
         sigc=sigi
      elseif(model.eq.7)then
#ifndef __QGSJETII__
         stop'please compile with requested model'
#endif
         NITER=5000
         if(idtarg.eq.0)then
            e0=dble(elab)
            icp=idtrafo('nxs','qgs',idproj)
            call qgini(e0,icp,maproj,matarg)
         endif
         CALL qgcrossc(NITER,GTOT,GPROD,GABS,GDD,GQEL,GCOH)
         sigt=sngl(GTOT)
         sigi=sngl(GPROD)
         sigc=sngl(GABS)
         sige=sigt-sigi
c     elseif(model.eq.9)then
c     sigt=0.
c     sige=0.
c     sigc=0.
c     sigi=flucrse(ekin,maproj,matarg,idtarg)
      else
         sigt=0.
         sigi=0.
         sige=0.
         sigc=0.
      endif
      end


      subroutine m2XXFZ( a,b)
      double precision a,b
#ifndef __QGSJET01__
      stop'please compile with requested model'
#endif
      CALL XXFZ(a,b)
      end

      subroutine m3SIGMA(ek,idpro,idtar,latar,matar,sigi,sige)
#ifndef __GHEISHA__
      stop'please compile with requested model'
#endif
      call ghecrse(ek,idpro,idtar,latar,matar,sigi,sige)
      end

      subroutine m6SIGMA(icl,engy,stot,sela,sine,sdifr,slela,Rho)
      dimension sdifr0(3)
#ifndef __SIBYLL__
      stop'please compile with requested model'
#endif
      if(icl.eq.1)then
         L=2
      elseif(icl.eq.2)then
         L=1
      else
         L=3
      endif
      call SIB_SIGMA_HP(L,engy,stot,sela,sine,sdifr0,slela,Rho)
      sdifr=sdifr0(1)+sdifr0(2)+sdifr0(3)
      end


      subroutine m7SIGMA(stot,scut,sine,slela)
      double precision GzZ0(5),pi,bm,amws
      common /qgarr1/  ia(2),icz,icp
      common /qgarr6/  pi,bm,amws
#ifndef __QGSJETII__
      stop'please compile with requested model'
#endif
      ia2save=ia(2)
      bmsave=bm
      ia(2)=1
      call qgfz(0.d0,gzz0,0,0)
      scut=sngl(gzz0(2))/2.     !cut pomerons cross-section
      stot=sngl(gzz0(1))        !tot cross-section
      sine=sngl(gzz0(2)+gzz0(3)+gzz0(4))/2. !inelastic cross section
      slela=sngl(gzz0(5))
      ia(2)=ia2save
      bm=bmsave
      end

      subroutine m8SIGMA(stot,scut,sine,sela,slela,ssd)
C     cross sections
      INTEGER IPFIL,IFAFIL,IFBFIL
      DOUBLE PRECISION SIGTOT,SIGELA,SIGVM,SIGINE,SIGNDF,SIGDIR,
     &     SIGLSD,SIGHSD,SIGLDD,SIGHDD,SIGCDF,
     &     SIGPOM,SIGREG,SIGHAR,SIGTR1,SIGTR2,SIGLOO,
     &     SIGDPO,SIG1SO,SIG1HA,SLOEL,SLOVM,SIGCOR,
     &     FSUP,FSUD,FSUH,ECMFIL,P2AFIL,P2BFIL
      COMMON /POCSEC/ SIGTOT,SIGELA,SIGVM(0:4,0:4),SIGINE,SIGNDF,SIGDIR,
     &     SIGLSD(2),SIGHSD(2),SIGLDD,SIGHDD,SIGCDF(0:4),
     &     SIGPOM,SIGREG,SIGHAR,SIGTR1(2),SIGTR2(2),SIGLOO,
     &     SIGDPO(4),SIG1SO,SIG1HA,SLOEL,SLOVM(4,4),SIGCOR,
     &     FSUP(2),FSUD(2),FSUH(2),ECMFIL,P2AFIL,P2BFIL,
     &     IPFIL,IFAFIL,IFBFIL
#ifndef __Phojet__
         stop'please compile with requested model'
#endif
      stot=sngl(SIGTOT)
      sine=sngl(SIGTOT-SIGELA)
      sela=sngl(SIGELA)
      slela=sngl(SLOEL)
      ssd=sngl(SIGLSD(1)+SIGHSD(1)+SIGLSD(2)+SIGHSD(2))
      sdd=sngl(SIGLDD+SIGHDD)
      scut=sine-ssd-sdd-sngl(SIGCDF(0))
      end

      subroutine m9SIGMA(stot,sine,sela)
c     include "epos.inc"
c     double precision PLA,EKIN1,SHPTOT,SHPINE,ZZ,AA,Sel,Zl
      stot=0.
      sine=0.
      sela=0.
c     if(iclpro.eq.1)then
c     IP=13
c     elseif(iclpro.eq.2)then
c     IP=1
c     else
c     IP=15
c     endif
c     PLA=dble(pnll)
c     EKIN1=dble(ekin)
c     ZZ=1
c     AA=1
c     stot=SHPTOT(IP,PLA)
c     sine=SHPINE(IP,PLA)
c     call SIGELH(IP,ZZ,AA,EKIN1,PLA,Sel,Zl)
c     sela=sngl(Sel)
      end

      subroutine decaymod(ip,iret)
      include 'epos.inc'
      if(model.eq.4)then
#ifndef __PYTHIA__
         stop'please compile with requested model'
#endif
        call emspyt(iret,ip)
      endif
      end
