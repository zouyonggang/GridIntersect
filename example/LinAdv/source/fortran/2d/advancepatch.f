c**********************************************************************
c 计算通量.
c**********************************************************************
      subroutine advancepatch2d(velocity, 
     &                          dt, 
     &                          num_edges,
     &                          num_extend_cells,
     &                          edge_coords, 
     &                          cell_extend_coords, 
     &                          eac_extent, 
     &                          eac_index, 
     &                          uval,
     &                          flux)
c**********************************************************************
      implicit none
c**********************************************************************
c input parameters:
      double precision velocity, dt

      ! 网格片内部边的数目.
      integer num_edges
      ! 网格片内部及第1层影像区中的单元数之和.
      integer num_extend_cells

      ! <边, 单元>拓扑.
      integer eac_extent(0:num_edges)
      integer eac_index(0:eac_extent(num_edges)-1)

      ! 边中心坐标.
      double precision edge_coords(0:num_edges*2-1)

      ! 单元中心坐标.
      double precision cell_extend_coords(0:num_extend_cells*2-1)

      ! 守恒量：中心量数据片
      double precision uval(0:num_extend_cells-1)

c output parameters:
      ! 通量：边心量数据片.
      double precision flux(0:num_edges-1)

c**********************************************************************
c     local parameters:
      integer i, idx1, idx2
      double precision tmp_coord
      double precision epsilon
      parameter (epsilon = 0.00000001)
c**********************************************************************

      do i = 0, num_edges - 1
        flux(i) = 0.0d0
      enddo

      do i = 0, num_edges - 1
        !用于识别的X, Y方向的边的辅助量.
        tmp_coord = edge_coords(2*i) - epsilon 

        ! 如果与一个边相邻的单元有两个, 则说明该边是内部边,
        if (eac_extent(i+1) - eac_extent(i) .eq. 2) then
          ! 取出与第i个边相邻的两个单元的索引.
          idx1 = eac_index(eac_extent(i))
          idx2 = eac_index(eac_extent(i)+1)

          ! 对于垂直于X坐标轴的边, 用其左边单元的物理解计算通量,
          ! 其它边上的通量值均为0. 
          if (tmp_coord .gt. cell_extend_coords(2*idx1)) then 
            flux(i) = dt * velocity * uval(idx1)
          else if (tmp_coord .gt. cell_extend_coords(2*idx2)) then   
            flux(i) = dt * velocity * uval(idx2)
          endif

        ! 如果与一个面相邻的单元只有一个，则说明该面是物理边界上的面.
        else if (eac_extent(i+1) - eac_extent(i) .eq. 1) then
          idx1 = eac_index(eac_extent(i))
          if (tmp_coord .gt. cell_extend_coords(2*idx1)) then
            ! X方向右边界上的边.
            flux(i) = dt * velocity * uval(idx1)
          endif
        else 
          write (6,*) "Error:wrong number of cells neighbouring to edge"
          stop
        endif
      enddo

      return
      end
