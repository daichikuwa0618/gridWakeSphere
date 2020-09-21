! ======================================================================
!   尾っぽをみたいグリッドをつくる
! ======================================================================

      program gridWake
      double precision, allocatable :: xa(:,:),ya(:,:) ! 格子点座標
      integer imax, jmax ! 格子点数
      double precision r, dr, theta, dTheta, rmin, drmin
      double precision pi, rHead, rTop, rTail, drFactor, dnstyRate
      double precision jGrowth, djGrouth, iGrowth, diGrowth, drTail

      ! 定数
      ! ====== input ======
      imax      = 501    ! 周方向 dnstyRateで割り切れる数+1 がいいかも
      jmax      = 501    ! 半径方向
      rmin      = 1.0d-2 ! 球の半径 [m]
      rHead     = rmin*2.d0 ! 球前方最大半径 [m]
      rTail     = rmin*16.d0 ! Wake 側の最大半径 [m]
      rTop      = rmin*4.d0 ! 頂点の半径
      dnstyRate = 1.0d0 ! 頭と後ろの密度比
      drFactor  = 400.d0 ! 最小格子幅の係数
      ! ====== input ここまで ======

      iTop      = nint(dble(imax)/(dnstyRate+1)) ! 頂点でのi
      drmin     = rmin/drFactor ! 最小格子幅
      pi        = dacos(-1.0d0)

      jGrowth = 0.d0 ! 半径方向の成長率
      djGrouth = 1.d-7 ! jGrowthの探索用
      iGrowth = 0.d0 ! 周方向の成長率
      djGrouth = 1.d-7 ! iGrowthの探索用

      allocate(xa(imax, jmax))
      allocate(ya(imax, jmax))

! jGrowth を求める
      ! ちょっと超えるまで
      do while (r.lt.rTop)
         jGrowth = jGrowth + djGrouth
         r = rmin
         do j=2,jmax
            dr = drmin*dexp(jGrowth*dble(j-2))
            r  = r + dr
         end do
      end do
      write(6,*)'1 jGrowth=',jGrowth
      djGrouth = djGrouth/1.0d2 ! djGrouth を小さくして精度を上げる
      ! ちょっと下回るまで
      do while (r.gt.rTop)
         jGrowth = jGrowth - djGrouth
         r = rmin
         do j=2,jmax
            dr = drmin*dexp(jGrowth*dble(j-2))
            r  = r + dr
         end do
      end do
      write(6,*)'2 jGrowth=',jGrowth

! ======================
!   領域全体に配置
! ======================
      r = rmin
      do j=1,jmax
         dr = drmin*dexp(jGrowth*dble(j-2))
         r  = r+dr
         ! 前方
         ! j につれて線形に楕円→円に近づく変数
         cirOvalHead = ((1-rTop/rHead)*j
     &             +(jmax*rTop/rHead-1))/(jmax-1)
         do i=1,iTop
            xa(i,j) = (rHead/rTop)*r*cirOvalHead
     &                 *dcos(((dble(i+1-2*iTop))/(dble(1-iTop)))*pi/2)
            ya(i,j) = r*dsin(((dble(i+1-2*iTop))/(dble(1-iTop)))*pi/2)
         end do
         ! 後方
         ! j につれて線形に円→楕円に近づく変数
         cirOvalTail = ((1-rTop/rTail)*j
     &             +(jmax*rTop/rTail-1))/(jmax-1)
         do i=iTop+1,imax
            xa(i,j) = (rTail/rTop)*r*cirOvalTail
     &                 *dcos(((dble(i-imax))/(dble(iTop-imax)))*pi/2)
            ya(i,j) = r*dsin(((dble(i-imax))/(dble(iTop-imax)))*pi/2)
         end do
      end do

! ===========================================
!   ファイルの作成
! ===========================================
      open(10, file='grid.dat', status='replace')
      write(10,*) imax,jmax
      do j=1,jmax
         do i=1,imax
            write(10,*) xa(i,j),ya(i,j)
         end do
      end do
      close(10)

 666  stop
      end program


