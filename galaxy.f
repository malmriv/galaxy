      program galaxy
        implicit none
        real*8, dimension (:,:), allocatable :: x,y,vx,vy,ax,ay,KE,UE
        real*8, dimension (:,:), allocatable :: vxcol, vycol
        real*8, dimension (:), allocatable :: m, R
        integer, dimension(:,:), allocatable :: vlist
        real*8 h, t, tmax, pos1, pos2, dp, dist, vradius
        real start, finish
        real*8 calc_x, calc_v, calc_a, aux, coll, total_KE, total_UE
        integer i, j, k, l, n_bodies, N, nonzero

        call cpu_time(start)
        n_bodies = 993
        N = 50
        tmax = 1.d0 !1 unit ~ current age of the sun
        h = tmax/float(N)
        vradius = 0.5 !In normalised distance units

c       Allocate memory for each vector
        allocate(x(1:n_bodies,1:N),y(1:n_bodies,1:N))
        allocate(vx(1:n_bodies,1:N),vy(1:n_bodies,1:N))
        allocate(vxcol(1:n_bodies,1:N),vycol(1:n_bodies,1:N))
        allocate(ax(1:n_bodies,1:N),ay(1:n_bodies,1:N))
        allocate(KE(1:n_bodies,1:N),UE(1:n_bodies,1:N))
        allocate(m(1:n_bodies),R(1:n_bodies))
        allocate(vlist(1:n_bodies,1:n_bodies))

c       Initialize the vectors.
        R = 0.d0
        UE = 0.d0
        KE = 0.d0
        t = 0.d0
        m = 0.d0
        x = 0.d0
        y = 0.d0
        vx = 0.d0
        vy = 0.d0
        ax = 0.d0
        ay = 0.d0
        dist = 0.d0

c       Data is generated and normalised by an R script.
        open(100,file="./dataset.txt",status='old')
        do i=1,n_bodies
            read(100,*) m(i), x(i,1), y(i,1), vx(i,1), vy(i,1), R(i)
        end do
        close(100)

c       We open different files for different magnitudes. This may seem
c       inefficient for a small number of bodies, but a simulation with
c       hundreds of bodies, each one evaluated thousands of times,
c       would quickly yield a very unmanageable .txt file if we were to
c       store everything in one go.
        call system("mkdir ./results")
        open(10,file="./results/x.txt")
        open(20,file="./results/y.txt")
        open(30,file="./results/vx.txt")
        open(40,file="./results/vy.txt")
        open(50,file="./results/ax.txt")
        open(60,file="./results/ay.txt")
        open(70,file="./results/t.txt")
        open(80,file="./results/energy.txt")

c       Evaluate accelerations to set initial conditions
        i = 1
        do j=1,n_bodies !j = body whose motion is being computed
          do k=1,n_bodies !k = its interactions with the k-th body
            if(j .ne. k) then !Avoid self-interactions
            ax(j,i) = ax(j,i) + calc_a(m(k),x(j,i),x(k,i),y(j,i),y(k,i))
            ay(j,i) = ay(j,i) + calc_a(m(k),y(j,i),y(k,i),x(j,i),x(k,i))
            end if
          end do
        end do

c       Compute initial Verlet list
        vlist = 0
        do j=1,n_bodies
        do k=1,n_bodies
          if(j.eq.k) cycle
          dist = (x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0
          if(dist .lt. vradius**2.0) vlist(j,k) = k
        end do
        call quicksort(vlist(j,1:n_bodies),1,n_bodies)
        call reverse(vlist(j,1:n_bodies),n_bodies)
        end do

c       Save the initial time
        write(70,*) t

c       Iterative process starts here!
        do i=2,N

        do j=1,n_bodies
c         Compute positions
          x(j,i) = calc_x(h,x(j,i-1),vx(j,i-1),ax(j,i-1))
          y(j,i) = calc_x(h,y(j,i-1),vy(j,i-1),ay(j,i-1))
        end do

c       Compute Verlet list ever 5 steps, best to determine later
          if(modulo(i,5) .eq. 0) then
            vlist = 0
            do j=1,n_bodies
            do k=1,n_bodies
              if(j.eq.k) cycle
              dist = (x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0
            if(dist .lt. vradius**2.0) then
              vlist(j,k) = k
            end if
            end do
            call quicksort(vlist(j,1:n_bodies),1,n_bodies)
            call reverse(vlist(j,1:n_bodies),n_bodies)
            end do
          end if

c         Compute accelerations
          do j=1,n_bodies
            if(nonzero(vlist(j,1:n_bodies),n_bodies) .eq. 0) cycle
            do k=1,nonzero(vlist(j,1:n_bodies),n_bodies)
              l = vlist(j,k)
              ax(j,i) = ax(j,i) + calc_a(m(l),x(j,i),x(l,i),
     &        y(j,i),y(l,i))
              ay(j,i) = ay(j,i) + calc_a(m(l),y(j,i),y(l,i),
     &        x(j,i),x(l,i))
            end do
          end do

          do j=1,n_bodies
c           Compute velocities
            vx(j,i) = calc_v(h,vx(j,i-1),ax(j,i-1),ax(j,i))
            vy(j,i) = calc_v(h,vy(j,i-1),ay(j,i-1),ay(j,i))
          end do

c         Save velocities to compute collisions
          vxcol = vx
          vycol = vy

          do j=1,(n_bodies-1)
            do k=(j+1),n_bodies
c           Compute velocities in case a collision has taken place.
c           Define the necessary parameters to check for collisions
              aux = dsqrt((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0)
              coll = dp(vx(k,i)-vx(j,i),vy(k,i)-vy(j,i),x(k,i)-x(j,i),
     &        y(k,i)-y(j,i))
              if((aux .le. (R(k)+R(j))) .and. (coll .le. 0.d0)) then
c             First body's velocity
                write(*,*) "Collision:",j,k
                vxcol(j,i) = vx(j,i) -2.*m(k)/(m(j)+m(k))
     &            * dp(vx(j,i)-vx(k,i),vy(j,i)-vy(k,i),
     &              x(j,i)-x(k,i),y(j,i)-y(k,i))
     &            / ((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0)
     &            * (x(j,i)-x(k,i))
                vycol(j,i) = vy(j,i) -2.*m(k)/(m(j)+m(k))
     &            * dp(vx(j,i)-vx(k,i),vy(j,i)-vy(k,i),
     &              x(j,i)-x(k,i),y(j,i)-y(k,i))
     &            / ((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0)
     &            * (y(j,i)-y(k,i))
c             Second body's velocity
                vxcol(k,i) = vx(k,i) -2.*m(j)/(m(j)+m(k))
     &            * (dp(vx(k,i)-vx(j,i),vy(k,i)-vy(j,i),
     &              x(k,i)-x(j,i),y(k,i)-y(j,i))
     &            / ((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0))
     &            * (x(k,i)-x(j,i))
                vycol(k,i) = vy(k,i) -2.*m(j)/(m(j)+m(k))
     &            *  dp(vx(k,i)-vx(j,i),vy(k,i)-vy(j,i),
     &               x(k,i)-x(j,i),y(k,i)-y(j,i))
     &            / ((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0)
     &            * (y(k,i)-y(j,i))
              end if
          end do
          end do

c         Set new velocities
          vx = vxcol
          vy = vycol

          do j = 1,n_bodies
c             Energies
              KE(j,i) = 0.5d0*m(j)*(vx(j,i)**2.d0+vy(j,i)**2.d0)
          end do


          do j=1,n_bodies
c           Compute the potential energy
            do k=j+1,n_bodies
              aux = dsqrt((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0)
              UE(j,i) = UE(j,i) - (m(j)*m(k))/aux
            end do
          end do


c         Update time
          t=t+h
          write(70,*) t
        end do

c       Save results. i-th row = i-th iteration, j-th column = j-th object
        do i=1,N
              write(10,*) (x(j,i),j=1,n_bodies)
              write(20,*) (y(j,i),j=1,n_bodies)
              write(30,*) (vx(j,i),j=1,n_bodies)
              write(40,*) (vy(j,i),j=1,n_bodies)
              write(50,*) (ax(j,i),j=1,n_bodies)
              write(60,*) (ay(j,i),j=1,n_bodies)
        end do

        do i=2,N
          total_KE = 0.d0
          total_UE = 0.d0
          do j=1,n_bodies
            total_KE = total_KE + KE(j,i)
            total_UE = total_UE + UE(j,i)
          end do
          write(80,*) i,total_KE,total_UE
        end do

c       Close every file
        do i=10,80,10
        close(i)
        end do
        call cpu_time(finish)
        write(*,*) "Total runtime (in seconds):",finish-start
        write(*,*) "Data saved in ./results/"
      end program galaxy

c     These functions evaluate angular position, velocity and accelera-
c     tion. Defining them here will make this program more adaptable.

      function calc_x(h,x,v,a) result(ans)
        implicit none
        real*8 h, x, v, a, ans
        ans = x + h*v + ((h**2.0)/2.0)*a
        return
      end function calc_x

      function calc_v(h,v,a1,a2) result(ans)
        implicit none
        real*8 h, v, a1, a2, ans
        ans = v + (h/2.0)*(a1+a2)
        return
      end function calc_v

      function calc_a(m,x1,x2,y1,y2) result(ans)
        implicit none
        real*8 m,x1,x2,y1,y2,dist,ans
        dist = dsqrt((x1-x2)**2.0+(y1-y2)**2.0)**3.0
        ans =  -m*(x1-x2)/dist
        return
      end function calc_a

c     Dot product function
      function dp(u1,u2,v1,v2) result(ans)
        implicit none
        real*8 u1,u2,v1,v2,ans
        ans = u1*v1+u2*v2
        return
      end function dp

c     Quicksort algorithm
      recursive subroutine quicksort(a, first, last)
        implicit none
        integer a(*), first, last, i, j, aux
        real*8 x, t

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
           do while (a(i) .lt. x)
              i=i+1
           end do
           do while (x .lt. a(j))
              j=j-1
           end do
           if (i .ge. j) exit
           t = a(i);  a(i) = a(j);  a(j) = t
           i=i+1
           j=j-1
        end do
        if (first < i-1) call quicksort(a, first, i-1)
        if (j+1 < last)  call quicksort(a, j+1, last)
      end subroutine quicksort

      subroutine reverse(a,length)
        implicit none
        integer a(*), i, length, aux
        do i=1,floor(length/2.0)
          aux = a(i)
          a(i) = a(length-i+1)
          a(length-i+1) = aux
        end do
      end subroutine reverse

      function nonzero(a,length) result(count)
        implicit none
        integer a(*), length, i, count
        count = 0
        do i=1,length
          if(a(i) .ne. 0) count = count+1
        end do
      end function nonzero
