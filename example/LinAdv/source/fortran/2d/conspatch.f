c**********************************************************************
c 更新守恒量.
c**********************************************************************
      subroutine conspatch2d(d_dx,
     &                       num_cells, 
     &                       num_edges,
     &                       cell_coords, 
     &                       edge_coords, 
     &                       cae_extent,
     &                       cae_index,
     &                       uval_current,
     &                       uval_new,
     &                       flux)
c**********************************************************************
      implicit none
c**********************************************************************
c input parameters:
      ! X方向网格步长. 
      double precision d_dx
      ! 网格片内部单元和边的数目.
      integer num_cells, num_edges

      ! 单元中心坐标.
      double precision cell_coords(0:num_cells*2-1)
      ! 边中心坐标.
      double precision edge_coords(0:num_edges*2-1)

      ! <单元, 边>相邻信息.
      integer cae_extent(0:num_cells)
      integer cae_index(0:cae_extent(num_cells)-1)

      ! 收恒量：中心量数据片
      double precision uval_current(0:num_cells-1)
      double precision uval_new(0:num_cells-1)
     
c output parameters:
      ! 通量：边心量数据片.
      double precision flux(0:num_edges-1)

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
        tmp_coord_1 = cell_coords(2*i) + epsilon;
        tmp_coord_2 = cell_coords(2*i) - epsilon;

        do j = cae_extent(i), cae_extent(i+1)-1
          if (tmp_coord_1 .lt. edge_coords(2*cae_index(j))) then
            ! 施加X方向右边的通量.
            tmp_uval = tmp_uval - flux(cae_index(j)) * recip_dx 
          endif

          if (tmp_coord_2 .gt. edge_coords(2*cae_index(j))) then
            ! 施加X方向左边的通量.
            tmp_uval = tmp_uval + flux(cae_index(j)) * recip_dx
            endif
         enddo

         uval_new(i) = tmp_uval
      enddo

      return
      end
               
