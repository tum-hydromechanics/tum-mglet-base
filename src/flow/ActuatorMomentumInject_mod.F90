MODULE ActuatorMomentumInject
   USE core_mod
   USE flowcore_mod
   IMPLICIT NONE(type, external)
   PRIVATE

   !public
   PUBLIC :: init_ActuatorMomentumInject, finish_ActuatorMomentumInject, calculate_ActuatorMomentumInject
   
   !  private
   REAL(realk), PRIVATE :: rCentre(2)           ! Coordinates of the Propeller rotation center (x, y)
   REAL(realk), PRIVATE :: Propeller_Z_coordination ! Propeller location in z-direction
   REAL(realk), PRIVATE :: rw                   ! Angular velocity (rad/s)
   LOGICAL, PROTECTED :: has_ActuatorMomentumInject = .FALSE.  !  ! activity flag 
   REAL(realk), PRIVATE :: inR                  ! Inner radius of the actuator zone
   REAL(realk), PRIVATE :: outR                 ! Outer radius of the actuator zone
   REAL(realk), PRIVATE :: area_dt              ! Time step used to form sector area (should be as same as rotate_dt)
   REAL(realk), PRIVATE :: rotate_dt            ! Time step used in rotation update
   REAL(realk), PRIVATE :: centrePoint(2)       ! Center point of the current blade element (x, y)
   REAL(realk), PRIVATE :: previous_startangle1 ! Start angle 1 from previous timestep (blade 1)
   REAL(realk), PRIVATE :: previous_startangle2 ! Start angle 2 from previous timestep (blade 2)
   REAL(realk), PRIVATE :: angle1               ! Current blade angle 1
   REAL(realk), PRIVATE :: angle2               ! Current blade angle 2
   REAL(realk), PRIVATE  :: momsrc_x, momsrc_y, momsrc_z   ! momentum source
   

contains
   SUBROUTINE init_ActuatorMomentumInject()
      WRITE(*,*) 'init_ActuatorMomentumInject called'
       ! leaving inactive if no parameters specified
    has_ActuatorMomentumInject = .FALSE.
    IF (.NOT. fort7%exists("/flow/ActuatorMomentumInject_info")) RETURN   
    ! retrieving parameters from parameters.json
    !centre position
    CALL fort7%get_array("/flow/ActuatorMomentumInject_info/rCentre", rCentre)
    !Z location
    CALL fort7%get_value("/flow/ActuatorMomentumInject_info/Propeller_Z_coordination", Propeller_Z_coordination)
    !angular velocity 
    CALL fort7%get_value("/flow/ActuatorMomentumInject_info/rw", rw)
    !inner radius 
    CALL fort7%get_value("/flow/ActuatorMomentumInject_info/inR", inR)
    !outer radius 
    CALL fort7%get_value("/flow/ActuatorMomentumInject_info/outR", outR)
    !area calculation time step 
    CALL fort7%get_value("/flow/ActuatorMomentumInject_info/area_dt", area_dt)
    !rotation time step 
    CALL fort7%get_value("/flow/ActuatorMomentumInject_info/rotate_dt", rotate_dt)
    !initial angles 
    CALL fort7%get_value("/flow/ActuatorMomentumInject_info/previous_startangle1", previous_startangle1)
    CALL fort7%get_value("/flow/ActuatorMomentumInject_info/previous_startangle2", previous_startangle2)
    ! momentum source
    CALL fort7%get_value("/flow/Momentum_Source/momsrc_x", momsrc_x)
    CALL fort7%get_value("/flow/Momentum_Source/momsrc_y", momsrc_y)
    CALL fort7%get_value("/flow/Momentum_Source/momsrc_z", momsrc_z)
      !calculate angles
       angle1 = previous_startangle1 + area_dt*rw
       angle2 = previous_startangle2 + area_dt*rw
      ! display obtained parameters
      IF (myid == 0) THEN
            WRITE(*, '("ActuatorMomentumInject_info TERM:")')
            WRITE(*, '(2X, "CentrePosition: ", 2(G0, 1X))') rCentre
            WRITE(*, '(2X, "AngularVelocity: ", 1(G0, 1X))') rw
            WRITE(*, '(2X, "InnerRadius: ", 1(G0, 1X))') inR
            WRITE(*, '(2X, "OuterRadius: ", 1(G0, 1X))') outR
            WRITE(*, '(2X, "AreaTimeStep: ", 1(G0, 1X))') area_dt
            WRITE(*, '(2X, "RotateTimeStep: ", 1(G0, 1X))') rotate_dt
            WRITE(*, '(2X, "InitialAngle1: ", 1(G0, 1X))') previous_startangle1
            WRITE(*, '(2X, "InitialAngle2: ", 1(G0, 1X))') previous_startangle2
            WRITE(*, '(2X, "Propeller_Z_coordination: ", 1(G0, 1X))') Propeller_Z_coordination
            WRITE(*, '(2X, "Momentum_source_x: ", 1(G0, 1X))') momsrc_x
            WRITE(*, '(2X, "Momentum_source_y: ", 1(G0, 1X))') momsrc_y
            WRITE(*, '(2X, "Momentum_source_z: ", 1(G0, 1X))') momsrc_z
            WRITE(*, '()')
        END IF

      ! set active
      has_ActuatorMomentumInject = .TRUE.

   END SUBROUTINE init_ActuatorMomentumInject

   SUBROUTINE finish_ActuatorMomentumInject

      ! revoking activity
      has_ActuatorMomentumInject = .FALSE.

      RETURN

   END SUBROUTINE finish_ActuatorMomentumInject


   !calculate trapezoid area
   FUNCTION trapezoid_area(points) result(area)
      implicit none
      real(realk), intent(in) :: points(4, 2)
      real(realk) :: area
      real(realk) :: sorted_points(4, 2)
      real(realk) :: top_points(2, 2), bottom_points(2, 2)
      real(realk) :: top_base, bottom_base, height
      integer(intk) :: i, j, max_idx, k
      real(realk) :: y_values(4)

      do i = 1, 4
         do j = 1, 2
            sorted_points(i, j) = points(i, j)
         end do
         y_values(i) = points(i, 2)
      end do

      do i = 1, 3
         max_idx = i
         do j = i + 1, 4
            if (sorted_points(j, 2) > sorted_points(max_idx, 2)) then
               max_idx = j
            end if
         end do
         if (max_idx /= i) then
            do k = 1, 2
               call swap(sorted_points(i, k), sorted_points(max_idx, k))
            end do
         end if
      end do

      top_points = sorted_points(1:2, :)
      bottom_points = sorted_points(3:4, :)

      top_base = abs(top_points(2, 1) - top_points(1, 1))
      bottom_base = abs(bottom_points(2, 1) - bottom_points(1, 1))

      height = (top_points(1, 2) + top_points(2, 2))/2.0_realk - (bottom_points(1, 2) + bottom_points(2, 2))/2.0_realk

      area = 0.5_realk*(top_base + bottom_base)*height

   contains
      SUBROUTINE swap(a, b)
         real(realk), intent(inout) :: a, b
         real(realk) :: temp
         temp = a
         a = b
         b = temp
      END SUBROUTINE swap
   END FUNCTION trapezoid_area

   !identify grid Vertices 
   FUNCTION find_adjacent(p1, p2) result(adjacent_points)
      implicit none
      integer(intk), intent(in) :: p1, p2
      integer(intk) :: grid(2, 2) = reshape([1, 2, 3, 4], [2, 2])
      integer(intk) :: adjacent_points(2)
      integer(intk) :: i, j, idx

      idx = 1

      do i = 1, 2
         do j = 1, 2
            if (grid(i, j) /= p1 .and. grid(i, j) /= p2) then
               adjacent_points(idx) = grid(i, j)
               idx = idx + 1
            end if
         end do
      end do

   END FUNCTION find_adjacent
   
   !calculate triangle area
   FUNCTION triangle_area(point1, point2, point3) result(areaa)
      implicit none
      real(realk), intent(in) :: point1(2), point2(2), point3(2)
      real(realk) :: areaa, base, height

      base = sqrt((point2(1) - point1(1))**2 + (point2(2) - point1(2))**2)

      height = sqrt((point2(1) - point3(1))**2 + (point2(2) - point3(2))**2)

      areaa = 0.5_realk*base*height
   END FUNCTION triangle_area

   !Detect intersection between a line segment and a circle
   FUNCTION line_circle_intersection(x1, y1, x2, y2, a, b, r) result(points)
      use, intrinsic :: ieee_arithmetic
      implicit none

      real(realk), intent(in) :: x1, y1, x2, y2, a, b, r
      real(realk), allocatable :: points(:, :)
      real(realk) :: dx, dy, coef_a, coef_b, coef_c, disc, sqrt_disc, t1, t2
      real(realk), allocatable :: temp_points(:, :)
      integer(intk) :: count
      real(realk), parameter :: tol = 1.0_realk-10

      dx = x2 - x1
      dy = y2 - y1

      coef_a = dx*dx + dy*dy

      
      if (coef_a < tol) then
         if (abs((x1 - a)**2 + (y1 - b)**2 - r**2) < tol) then
            allocate (points(1, 2))
            points(1, 1) = x1
            points(1, 2) = y1
         else
            allocate (points(1, 2))
            points = reshape([ieee_value(0.0_realk, ieee_quiet_nan), ieee_value(0.0_realk, ieee_quiet_nan)], [1, 2])

         end if
         return
      end if

      coef_b = 2.0_realk*(dx*(x1 - a) + dy*(y1 - b))
      coef_c = (x1 - a)**2 + (y1 - b)**2 - r**2
      disc = coef_b**2 - 4.0_realk*coef_a*coef_c

      if (disc < -tol) then
         allocate (points(1, 2))
         points = reshape([ieee_value(0.0_realk, ieee_quiet_nan), ieee_value(0.0_realk, ieee_quiet_nan)], [1, 2])

         return
      end if

      allocate (temp_points(2, 2))
      count = 0

      if (disc >= -tol) then
         sqrt_disc = sqrt(max(disc, 0.0_realk))
         t1 = (-coef_b - sqrt_disc)/(2.0_realk*coef_a)
         t2 = (-coef_b + sqrt_disc)/(2.0_realk*coef_a)

         if (t1 >= -tol .and. t1 <= 1.0_realk + tol) then
            count = count + 1
            temp_points(count, 1) = x1 + t1*dx
            temp_points(count, 2) = y1 + t1*dy
         end if
         if (abs(t2 - t1) > tol .and. t2 >= -tol .and. t2 <= 1.0_realk + tol) then
            count = count + 1
            temp_points(count, 1) = x1 + t2*dx
            temp_points(count, 2) = y1 + t2*dy
         end if
      end if

      if (count == 0) then
         allocate (points(1, 2))
         points = reshape([ieee_value(0.0_realk, ieee_quiet_nan), ieee_value(0.0_realk, ieee_quiet_nan)], [1, 2])

      else
         allocate (points(count, 2))
         points = temp_points(1:count, :)
      end if

   END FUNCTION line_circle_intersection
   
   !Detect intersection between two line segment 
   FUNCTION lineSegmentIntersection(x1, y1, x2, y2, x3, y3, x4, y4) result(intersection)
      use, intrinsic :: ieee_arithmetic
      implicit none
      real(realk), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
      real(realk) :: dx1, dy1, dx2, dy2, denom, t, u, px, py
      real(realk), allocatable :: intersection(:)

      dx1 = x2 - x1
      dy1 = y2 - y1
      dx2 = x4 - x3
      dy2 = y4 - y3

      denom = dx1*dy2 - dy1*dx2

      if (abs(denom) > 1e-10) then  ! Avoid division by zero
         t = ((x3 - x1)*dy2 - (y3 - y1)*dx2)/denom
         u = ((x3 - x1)*dy1 - (y3 - y1)*dx1)/denom

         if (t >= 0.0_realk .and. t <= 1.0_realk .and. u >= 0.0_realk .and. u <= 1.0_realk) then
            px = x1 + t*dx1
            py = y1 + t*dy1
            allocate (intersection(2))
            intersection(1) = px
            intersection(2) = py
         else
            allocate (intersection(2))
            intersection = [ieee_value(0.0_realk, ieee_quiet_nan), ieee_value(0.0_realk, ieee_quiet_nan)]
         end if
      else
         allocate (intersection(2))
         intersection = [ieee_value(0.0_realk, ieee_quiet_nan), ieee_value(0.0_realk, ieee_quiet_nan)]
      end if
   END FUNCTION lineSegmentIntersection
   
   !Detect vertices in sector 
   LOGICAL FUNCTION PointInSector(px, py, centre, innerR, outerR, startAngle, endAngle)
      implicit none
      real(realk), intent(in) :: px, py
      real(realk), intent(in) :: centre(2)
      real(realk), intent(in) :: innerR, outerR, startAngle, endAngle
      real(realk) :: distance, angle
      real(realk) :: modStartAngle, modEndAngle

      modStartAngle = mod(startAngle, 2.0_realk*3.141592653589793_realk)
      modEndAngle = mod(endAngle, 2.0_realk*3.141592653589793_realk)

      distance = sqrt((px - centre(1))**2 + (py - centre(2))**2)

      angle = atan2(py - centre(2), px - centre(1))

      if (angle < 0.0_realk) then
         angle = angle + 2.0_realk*3.141592653589793_realk
      end if

      if (distance >= innerR .and. distance <= outerR .and. &
          ((modStartAngle <= modEndAngle .and. angle >= modStartAngle .and. angle <= modEndAngle) .or. &
           (modStartAngle > modEndAngle .and. (angle >= modStartAngle .or. angle <= modEndAngle)))) then
         PointInSector = .true.
      else
         PointInSector = .false.
      end if

   END FUNCTION PointInSector
   
   !calculate velocity after mommentum injection 
   SUBROUTINE calculate_ActuatorMomentumInject(ittot, uo, vo, wo)
      INTEGER(intk), INTENT(IN) :: ittot
      TYPE(field_t), INTENT(INOUT) :: uo, vo, wo
      INTEGER(intk) :: i, igrid, kk, jj, ii
      REAL(realk), ALLOCATABLE :: momentum_inject_area_3d(:, :, :) ! 3D array storing swept area fraction
      REAL(realk), POINTER, CONTIGUOUS :: uo_arr(:,:,:), vo_arr(:,:,:), wo_arr(:,:,:)
      
      DO i = 1, nmygrids
          igrid = mygrids(i)
          CALL get_mgdims(kk, jj, ii, igrid)
           IF (ALLOCATED(momentum_inject_area_3d)) DEALLOCATE(momentum_inject_area_3d)
           ALLOCATE(momentum_inject_area_3d(kk, jj, ii))
           momentum_inject_area_3d = 0.0_realk
           
          CALL uo%get_ptr(uo_arr, igrid)
          CALL vo%get_ptr(vo_arr, igrid)
          CALL wo%get_ptr(wo_arr, igrid)
          CALL Area_Calculation_single_grid(ittot, momentum_inject_area_3d, igrid)
          CALL inject_momentum_source(kk, jj, ii, uo_arr, vo_arr, wo_arr, momentum_inject_area_3d)
          
          IF (ALLOCATED(momentum_inject_area_3d)) DEALLOCATE(momentum_inject_area_3d)
      END DO
   END SUBROUTINE calculate_ActuatorMomentumInject
   
   !calculate area in every grid,output 3D matrix
   SUBROUTINE Area_Calculation_single_grid(step, momentum_inject_area_3d, igrid)
      use, intrinsic :: ieee_arithmetic
      INTEGER(intk), INTENT(IN) :: step
      REAL(realk), INTENT(INOUT), ALLOCATABLE :: momentum_inject_area_3d(:, :, :)
      INTEGER(intk), INTENT(IN) :: igrid
      
      
      real(realk)::sector1_startAngle, sector1_endAngle
      real(realk)::sector2_startAngle, sector2_endAngle
      real(realk)::add_angle, t
      REAL(realk) :: innerStartX1, innerStartY1, innerEndX1, innerEndY1
      REAL(realk) :: OuterStartX1, OuterStartY1, OuterEndX1, OuterEndY1
      REAL(realk) :: innerStartX2, innerStartY2, innerEndX2, innerEndY2
      REAL(realk) :: OuterStartX2, OuterStartY2, OuterEndX2, OuterEndY2
      INTEGER(intk) :: j, k
      INTEGER(intk) :: ni, nj, kk, jj, ii
      real(realk) :: x1, y1, x2, y2, x3, y3, x4, y4
      real(realk) :: xpoints(4), ypoints(4)
      integer(intk) :: q(4)
      integer(intk) :: mapping(4, 2)
      integer(intk) :: m, idx, New_m, m_count, New_idx(2)
      integer(intk) :: adjacent_points(2)
      real(realk) :: orthoPoint(2)
      integer(intk) :: calcupoint(2)
      integer(intk) :: New_New_m, New_New_idx
      real(realk) :: New_orthoPoint(2)
      integer(intk) :: New_calcupoint(2)
      !A,B,C,D,E,F are intersection points
      real(realk) :: coordinateA(2), coordinateB(2), coordinateC(2), coordinateD(2)
      real(realk) :: coordinateE(2), coordinateF(2)
      real(realk) :: coordinate1(2), coordinate2(2), coordinate3(2), coordinate4(2)
      real(realk) :: coordinate7(2), coordinate8(2), coordinate9(2), coordinate10(2)
      real(realk) :: coordinate13(2), coordinate14(2), coordinate15(2), coordinate16(2)
      real(realk) :: coordinate19(2), coordinate20(2), coordinate21(2), coordinate22(2)
      real(realk) :: coordinate25(2), coordinate26(2), coordinate27(2), coordinate28(2)
      real(realk) :: coordinate31(2), coordinate32(2), coordinate33(2), coordinate34(2)
      real(realk), allocatable :: coordinate5(:, :), coordinate6(:, :)
      real(realk), allocatable :: coordinate11(:, :), coordinate12(:, :)
      real(realk), allocatable :: coordinate17(:, :), coordinate18(:, :)
      real(realk), allocatable :: coordinate23(:, :), coordinate24(:, :)
      real(realk), allocatable :: coordinate29(:, :), coordinate30(:, :)
      real(realk), allocatable :: coordinate35(:, :), coordinate36(:, :)
      real(realk) ::area, New_area
      real(realk), dimension(4, 2) :: trapezoid_Point
      real(realk) :: trapezoid_Area_calculated
      
      REAL(realk), ALLOCATABLE :: Xgrid(:, :), Ygrid(:, :)
      REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
      REAL(realk) :: dx, dy
      
      TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
      REAL(realk), POINTER, CONTIGUOUS :: dx_ptr(:), dy_ptr(:), dz_ptr(:)
      REAL(realk), POINTER, CONTIGUOUS :: ddx_ptr(:), ddy_ptr(:), ddz_ptr(:)
      REAL(realk), ALLOCATABLE :: momentum_inject_area(:, :)  ! 2D array storing swept area fraction 
      INTEGER(intk) :: z_layer
   
      !check status
      IF (.NOT. has_ActuatorMomentumInject) RETURN
   
      t = step*rotate_dt
      add_angle = rw*t
      centrePoint = rCentre
   
      !sectror1
      sector1_startAngle = previous_startangle1 + add_angle
      sector1_endAngle = angle1 + add_angle
      !sector2
      sector2_startAngle = previous_startangle2 + add_angle
      sector2_endAngle = angle2 + add_angle
      
      !endPoints-sector1
      innerStartX1 = centrePoint(1) + inR*cos(sector1_startAngle)
      innerStartY1 = centrePoint(2) + inR*sin(sector1_startAngle)
      innerEndX1 = centrePoint(1) + inR*cos(sector1_endAngle)
      innerEndY1 = centrePoint(2) + inR*sin(sector1_endAngle)
      OuterStartX1 = centrePoint(1) + outR*cos(sector1_startAngle)
      OuterStartY1 = centrePoint(2) + outR*sin(sector1_startAngle)
      OuterEndX1 = centrePoint(1) + outR*cos(sector1_endAngle)
      OuterEndY1 = centrePoint(2) + outR*sin(sector1_endAngle)
      
      !endPoints-sector2
      innerStartX2 = centrePoint(1) + inR*cos(sector2_startAngle)
      innerStartY2 = centrePoint(2) + inR*sin(sector2_startAngle)
      innerEndX2 = centrePoint(1) + inR*cos(sector2_endAngle)
      innerEndY2 = centrePoint(2) + inR*sin(sector2_endAngle)
      OuterStartX2 = centrePoint(1) + outR*cos(sector2_startAngle)
      OuterStartY2 = centrePoint(2) + outR*sin(sector2_startAngle)
      OuterEndX2 = centrePoint(1) + outR*cos(sector2_endAngle)
      OuterEndY2 = centrePoint(2) + outR*sin(sector2_endAngle)
      
      !form grid matrix
      CALL get_mgdims(kk, jj, ii, igrid)
      ni = ii
      nj = jj
      
      CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
      
      CALL get_field(dx_f, "DX")
      CALL get_field(dy_f, "DY")
      CALL get_field(dz_f, "DZ")
      CALL get_field(ddx_f, "DDX")
      CALL get_field(ddy_f, "DDY")
      CALL get_field(ddz_f, "DDZ")
      
      CALL dx_f%get_ptr(dx_ptr, igrid)
      CALL dy_f%get_ptr(dy_ptr, igrid)
      CALL dz_f%get_ptr(dz_ptr, igrid)
      CALL ddx_f%get_ptr(ddx_ptr, igrid)
      CALL ddy_f%get_ptr(ddy_ptr, igrid)
      CALL ddz_f%get_ptr(ddz_ptr, igrid)
      
      ALLOCATE(Xgrid(nj, ni), Ygrid(nj, ni))
      
      DO j = 1, nj
         DO k = 1, ni
            Xgrid(j, k) = minx + REAL(k - 1, realk) * dx_ptr(k)
            Ygrid(j, k) = miny + REAL(j - 1, realk) * dy_ptr(j)
         END DO
      END DO
      
      IF (.NOT. ALLOCATED(momentum_inject_area)) THEN
         ALLOCATE(momentum_inject_area(nj-1, ni-1))
      END IF
      momentum_inject_area = 0.0_realk
      
      
      DO j = 3, nj-2
         DO k = 3, ni-2
            x1 = Xgrid(j, k); y1 = Ygrid(j, k)
            x2 = Xgrid(j, k + 1); y2 = Ygrid(j, k + 1)
            x3 = Xgrid(j + 1, k); y3 = Ygrid(j + 1, k)
            x4 = Xgrid(j + 1, k + 1); y4 = Ygrid(j + 1, k + 1)
            
            xpoints = [x1, x2, x3, x4]
            ypoints = [y1, y2, y3, y4]
            
            q = [0, 0, 0, 0]
            mapping(1, :) = [2, 3]
            mapping(2, :) = [1, 4]
            mapping(3, :) = [1, 4]
            mapping(4, :) = [2, 3]
            
            Do m = 1, 4
               if (PointInSector(xpoints(m), ypoints(m), centrePoint, inR, outR, sector1_startAngle, sector1_endAngle) .or. &
                   PointInSector(xpoints(m), ypoints(m), centrePoint, inR, outR, sector2_startAngle, sector2_endAngle)) then
                  q(m) = 1
               else
                  q(m) = 0
               end if
            end do
            
            
            ! 4 vertices all inside sector
            if (sum(q) == 4) then
               momentum_inject_area(j, k) = 1

               ! 3 vertices inside sector
            elseif (sum(q) == 3) then
               
               do m = 1, 4
                     if (q(m) == 0) then
                        idx = m
                        exit
                     end if
                  end do
                  orthoPoint(1) = xpoints(idx)
                  orthoPoint(2) = ypoints(idx)
                  calcupoint = mapping(idx, :)
                  !Calculate A Point
                  coordinate1=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
                  if (.not. ieee_is_nan(coordinate1(1))) then
                     coordinateA = coordinate1
                  else
                     coordinate2=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                     if (.not. ieee_is_nan(coordinate2(1))) then
                        coordinateA = coordinate2
                     else
                        coordinate3=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                        if (.not. ieee_is_nan(coordinate3(1))) then
                           coordinateA = coordinate3
                        else
                         coordinate4=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                           if (.not. ieee_is_nan(coordinate4(1))) then
                              coordinateA = coordinate4
                           else
                              coordinate5 = line_circle_intersection(xpoints(calcupoint(1)), ypoints(calcupoint(1)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), inR)
                              if (.not. ieee_is_nan(coordinate5(1, 1))) then
                                 coordinateA = coordinate5(1, :)
                              else
                                 coordinate6 = line_circle_intersection(xpoints(calcupoint(1)), ypoints(calcupoint(1)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), outR)
                                 if (.not. ieee_is_nan(coordinate6(1, 1))) then
                                    coordinateA = coordinate6(1, :)
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
                  !Calculate B Point
                  coordinate7=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
                  if (.not. ieee_is_nan(coordinate7(1))) then
                     coordinateB = coordinate7
                  else
                     coordinate8=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                     if (.not. ieee_is_nan(coordinate8(1))) then
                        coordinateB = coordinate8
                     else
                     coordinate9=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                        if (.not. ieee_is_nan(coordinate9(1))) then
                           coordinateB = coordinate9
                        else
                     coordinate10=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                           if (.not. ieee_is_nan(coordinate10(1))) then
                              coordinateB = coordinate10
                           else
                              coordinate11 = line_circle_intersection(xpoints(calcupoint(2)), ypoints(calcupoint(2)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), inR)
                              if (.not. ieee_is_nan(coordinate11(1, 1))) then
                                 coordinateB = coordinate11(1, :)
                              else
                                 coordinate12 = line_circle_intersection(xpoints(calcupoint(2)), ypoints(calcupoint(2)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), outR)
                                 if (.not. ieee_is_nan(coordinate12(1, 1))) then
                                    coordinateB = coordinate12(1, :)
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
                  !calculate cut area
                  area = triangle_area(coordinateA, orthoPoint, coordinateB)
                  momentum_inject_area(j, k) = ((dx_ptr(k)*dy_ptr(j)) - area)/(dx_ptr(k)*dy_ptr(j))
                       ! Limit the result to reasonable range
                  IF (momentum_inject_area(j, k) < 0.0_realk) momentum_inject_area(j, k) = 0.0_realk
                  IF (momentum_inject_area(j, k) > 1.0_realk) momentum_inject_area(j, k) = 1.0_realk

                  ! 2 vertices inside sector
            elseif (sum(q) == 2) then
               
               m_count = 0
                  do New_m = 1, 4
                     if (q(New_m) == 0) then
                        m_count = m_count + 1
                        New_idx(m_count) = New_m
                        if (m_count == 2) exit
                     end if
                  end do
                  adjacent_points = find_adjacent(New_idx(1), New_idx(2))
                  !Calculate C Point
                  coordinate13=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
                  if (.not. ieee_is_nan(coordinate13(1))) then
                     coordinateC = coordinate13
                  else
                     coordinate14=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                     if (.not. ieee_is_nan(coordinate14(1))) then
                        coordinateC = coordinate14
                     else
                     coordinate15=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                        if (.not. ieee_is_nan(coordinate15(1))) then
                           coordinateC = coordinate15
                        else
                     coordinate16=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                           if (.not. ieee_is_nan(coordinate16(1))) then
                              coordinateC = coordinate16
                           else
                         coordinate17=line_circle_intersection(xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),xpoints(New_idx(1)),ypoints(New_idx(1)),centrePoint(1), centrePoint(2), inR)
                              if (.not. ieee_is_nan(coordinate17(1, 1))) then
                                 coordinateC = coordinate17(1, :)
                              else
                             coordinate18=line_circle_intersection(xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),xpoints(New_idx(1)),ypoints(New_idx(1)),centrePoint(1), centrePoint(2), outR)
                                 if (.not. ieee_is_nan(coordinate18(1, 1))) then
                                    coordinateC = coordinate18(1, :)
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
                  !Calculate D Point
                  coordinate19=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
                  if (.not. ieee_is_nan(coordinate19(1))) then
                     coordinateD = coordinate19
                  else
                     coordinate20=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                     if (.not. ieee_is_nan(coordinate20(1))) then
                        coordinateD = coordinate20
                     else
                     coordinate21=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                        if (.not. ieee_is_nan(coordinate21(1))) then
                           coordinateD = coordinate21
                        else
                     coordinate22=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                           if (.not. ieee_is_nan(coordinate22(1))) then
                              coordinateD = coordinate22
                           else
                         coordinate23=line_circle_intersection(xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),xpoints(New_idx(2)),ypoints(New_idx(2)),centrePoint(1), centrePoint(2), inR)
                              if (.not. ieee_is_nan(coordinate23(1, 1))) then
                                 coordinateD = coordinate23(1, :)
                              else
                             coordinate24=line_circle_intersection(xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),xpoints(New_idx(2)),ypoints(New_idx(2)),centrePoint(1), centrePoint(2), outR)
                                 if (.not. ieee_is_nan(coordinate24(1, 1))) then
                                    coordinateD = coordinate24(1, :)
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
                  !Calculate trapzoid area
                  trapezoid_Point(1, 1) = xpoints(New_idx(1))
                  trapezoid_Point(1, 2) = ypoints(New_idx(1))
                  trapezoid_Point(2, 1) = xpoints(New_idx(2))
                  trapezoid_Point(2, 2) = ypoints(New_idx(2))
                  trapezoid_Point(3, 1) = coordinateC(1)
                  trapezoid_Point(3, 2) = coordinateC(2)
                  trapezoid_Point(4, 1) = coordinateD(1)
                  trapezoid_Point(4, 2) = coordinateD(2)
                  trapezoid_Area_calculated = trapezoid_area(trapezoid_Point)
                  momentum_inject_area(j, k) = ((dx_ptr(k)*dy_ptr(j)) - trapezoid_Area_calculated)/(dx_ptr(k)*dy_ptr(j))
                       ! Limit the result to reasonable range
                  IF (momentum_inject_area(j, k) < 0.0_realk) momentum_inject_area(j, k) = 0.0_realk
                  IF (momentum_inject_area(j, k) > 1.0_realk) momentum_inject_area(j, k) = 1.0_realk

                  ! 1 vertices inside sector
            elseif (sum(q) == 1) then
               
                do New_New_m = 1, 4
                     if (q(New_New_m) == 1) then
                        New_New_idx = New_New_m
                        exit
                     end if
                  end do
                  New_orthoPoint(1) = xpoints(New_New_idx)
                  New_orthoPoint(2) = ypoints(New_New_idx)
                  New_calcupoint = mapping(New_New_idx, :)
                  !Calculate E point
                  coordinate25=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
                  if (.not. ieee_is_nan(coordinate25(1))) then
                     coordinateE = coordinate25
                  else
                      coordinate26=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                     if (.not. ieee_is_nan(coordinate26(1))) then
                        coordinateE = coordinate26
                     else
                     coordinate27=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                        if (.not. ieee_is_nan(coordinate27(1))) then
                           coordinateE = coordinate27
                        else
                     coordinate28=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                           if (.not. ieee_is_nan(coordinate28(1))) then
                              coordinateE = coordinate28
                           else
                    coordinate29=line_circle_intersection(xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),xpoints(New_New_idx),ypoints(New_New_idx), centrePoint(1), centrePoint(2), inR) 
                              if (.not. ieee_is_nan(coordinate29(1, 1))) then
                                 coordinateE = coordinate29(1, :)
                              else
                                  coordinate30=line_circle_intersection(xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)), xpoints(New_New_idx),ypoints(New_New_idx),centrePoint(1), centrePoint(2), outR)
                                 if (.not. ieee_is_nan(coordinate30(1, 1))) then
                                    coordinateE = coordinate30(1, :)
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
                  !Calculate F point
                  coordinate31=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
               if (.not. ieee_is_nan(coordinate31(1))) then
                  coordinateF = coordinate31
               else
                   coordinate32=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                  if (.not. ieee_is_nan(coordinate32(1))) then
                     coordinateF = coordinate32
                  else
                  coordinate33=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                     if (.not. ieee_is_nan(coordinate33(1))) then
                        coordinateF = coordinate33
                     else
                  coordinate34=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                        if (.not. ieee_is_nan(coordinate34(1))) then
                           coordinateF = coordinate34
                        else
                 coordinate35=line_circle_intersection(xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),xpoints(New_New_idx),ypoints(New_New_idx), centrePoint(1), centrePoint(2), inR) 
                           if (.not. ieee_is_nan(coordinate35(1, 1))) then
                              coordinateF = coordinate35(1, :)
                           else
                               coordinate36=line_circle_intersection(xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)), xpoints(New_New_idx),ypoints(New_New_idx),centrePoint(1), centrePoint(2), outR)
                              if (.not. ieee_is_nan(coordinate36(1, 1))) then
                                 coordinateF = coordinate36(1, :)
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
                  New_area = triangle_area(coordinateE, New_orthoPoint, coordinateF)
                  momentum_inject_area(j, k) = (New_area)/(dx_ptr(k)*dy_ptr(j))
                   ! Limit the result to reasonable range
                  IF (momentum_inject_area(j, k) < 0.0_realk) momentum_inject_area(j, k) = 0.0_realk
                  IF (momentum_inject_area(j, k) > 1.0_realk) momentum_inject_area(j, k) = 1.0_realk
            end if
         END DO
      END DO
      
      DEALLOCATE(Xgrid, Ygrid)
      
      ! form 3D matrix
      z_layer = 0
      DO k = 1, kk
         IF (minz + (dz_ptr(k) * (k - 1)) <= Propeller_Z_coordination .AND. &
             minz + (dz_ptr(k) * k) > Propeller_Z_coordination) THEN
            z_layer = k
            EXIT
         END IF
      END DO
      
      
      
      
      DO j = 3, nj-2
         DO k = 3, ni-2
            momentum_inject_area_3d(z_layer, j, k) = momentum_inject_area(j, k)
         END DO
      END DO
   
    IF (ALLOCATED(momentum_inject_area)) DEALLOCATE(momentum_inject_area)

   END SUBROUTINE Area_Calculation_single_grid
   
   !inject momentum in every grid
 SUBROUTINE inject_momentum_source(kk, jj, ii, uo, vo, wo, momentum_inject_area_3d)     
      INTEGER, INTENT(IN) :: kk, jj, ii
      REAL(realk), INTENT(INOUT) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)
      
      REAL(realk), INTENT(IN) :: momentum_inject_area_3d(kk, jj, ii)
      INTEGER :: k, j, i
      REAL(realk) :: factor_x, factor_y, factor_z
      
         IF (ANY(ISNAN(momentum_inject_area_3d))) THEN
         PRINT *, "ERROR: momentum_inject_area_3d contains NaNs"
         STOP
         END IF
       !debug information
       !WRITE(*,*) 'inject: kk=', kk, 'jj=', jj, 'ii=', ii
       !WRITE(*,*) 'inject: momsrc_x=', momsrc_x, 'momsrc_y=', momsrc_y, 'momsrc_z=', momsrc_z
       !WRITE(*,*) 'inject: momentum_inject_area_3d max=', MAXVAL(momentum_inject_area_3d), 'min=', MINVAL(momentum_inject_area_3d)
       !WRITE(*,*) 'inject: uo max=', MAXVAL(uo), 'min=', MINVAL(uo)
       !WRITE(*,*) 'inject: vo max=', MAXVAL(vo), 'min=', MINVAL(vo)
       !WRITE(*,*) 'inject: wo max=', MAXVAL(wo), 'min=', MINVAL(wo)
     DO i = 3, ii-2
        DO j = 3, jj-2
           DO k = 3, kk-2
              factor_x = 0.5 * (momentum_inject_area_3d(k, j, i) + momentum_inject_area_3d(k, j, i+1))
              factor_y = 0.5 * (momentum_inject_area_3d(k, j, i) + momentum_inject_area_3d(k, j+1, i))
              factor_z = 0.5 * (momentum_inject_area_3d(k, j, i) + momentum_inject_area_3d(k+1, j, i))
               
              uo(k, j, i) = uo(k, j, i) + (momsrc_x * factor_x) / rho
              vo(k, j, i) = vo(k, j, i) + (momsrc_y * factor_y) / rho
              wo(k, j, i) = wo(k, j, i) + (momsrc_z * factor_z) / rho
            END DO
        END DO
      END DO
      !WRITE(*,*) 'inject after: uo max=', MAXVAL(uo), 'min=', MINVAL(uo)
      !WRITE(*,*) 'inject after: vo max=', MAXVAL(vo), 'min=', MINVAL(vo)
      !WRITE(*,*) 'inject after: wo max=', MAXVAL(wo), 'min=', MINVAL(wo)
END SUBROUTINE inject_momentum_source



END MODULE 	ActuatorMomentumInject 