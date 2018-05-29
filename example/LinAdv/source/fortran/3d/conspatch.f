c**********************************************************************
c 更新守恒量.
c**********************************************************************
      subroutine conspatch3d(d_dx,
     &                       num_cells, 
     &                       num_faces,
     &                       cell_coords, 
     &                       face_coords, 
     &                       caf_extent,
     &                       caf_index,
     &                       uval_current,
     &                       uval_new,
     &                       flux)
c**********************************************************************
      implicit none
c**********************************************************************
c input parameters:
      ! X方向网格步长. 
      double precision d_dx
      ! 本地网格单元和面的数目.
      integer num_cells, num_faces

      ! 单元中心坐标.
      double precision cell_coords(0:num_cells*3-1)
      ! 边中心坐标.
      double precision face_coords(0:num_faces*3-1)

      ! <单元, 边>相邻信息.
      integer caf_extent(0:num_cells)
      integer caf_index(0:caf_extent(num_cells)-1)

      ! 守恒量：中心量数据片
      double precision uval_current(0:num_cells-1)
      double precision uval_new(0:num_cells-1)

c output parameters:
      ! 通量：面心量数据片.
      double precision flux(0:num_faces-1)

c**********************************************************************
c     local parameters:
      integer i, j
      double precision tmp_uval, tmp_coord_1, tmp_coord_2, recip_dx
      double precision epsilon
      parameter (epsilon = 0.00000001)
c**********************************************************************

      recip_dx = 1 / d_dx

      do i = 0, num_cells-1
        tmp_uval    = uval_current(i)
        tmp_coord_1 = cell_coords(3*i) + epsilon
        tmp_coord_2 = cell_coords(3*i) - epsilon

        do j = caf_extent(i), caf_extent(i+1)-1
          if (tmp_coord_1 .lt. face_coords(3*caf_index(j))) then
            ! 施加X方向右边的通量.
            tmp_uval = tmp_uval - flux(caf_index(j)) * recip_dx
          endif

          if (tmp_coord_2 .gt. face_coords(3*caf_index(j))) then
            ! 施加X方向右边的通量.
            tmp_uval = tmp_uval + flux(caf_index(j)) * recip_dx
          endif
        enddo

        uval_new(i) = tmp_uval
      enddo

      return
      end
               
