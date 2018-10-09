c**********************************************************************
c 初始化uval网格片数据.
c**********************************************************************
      subroutine initsetvalue(ndim,
     &                        d_dx,
     &                        num_cells,
     &                        num_nodes,
     &                        cell_coords,
     &                        node_coords,
     &                        can_extent,
     &                        can_index,
     &                        uval_current,
     &                        label_current)
c**********************************************************************
      implicit none
c**********************************************************************
c input parameters:
      double precision pi
      parameter(pi=3.1415926)

      ! 问题维数.
      integer ndim
      ! 网格片内部网格单元数.
      integer num_cells
      ! 网格片内部结点数.
      integer num_nodes

      ! 网格片内部网格单元中心坐标.
      double precision cell_coords(0:num_cells*ndim-1)
      ! 网格片内部结点坐标.
      double precision node_coords(0:num_nodes*ndim-1)

      !<单元, 结点>拓扑信息.
      integer can_extent(0:num_cells)
      integer can_index(0:can_extent(num_cells)-1)

c output parameters:
      ! X方向网格步长.
      double precision d_dx
      ! 中心量数据片.
      double precision uval_current(0:num_cells-1)
      double precision label_current(0:num_cells-1)

c**********************************************************************
c     local parameters:
      integer i, j, k
c**********************************************************************

      do i = 0, num_cells-1
        if (cell_coords(ndim*i) .le.
     &       (2.0d0-dsin(pi*cell_coords(ndim*i+1)))) then
          uval_current(i) = 2.0d0
        else
          uval_current(i) = 0.3d0
        endif
        !write (6,*) "cell_coords:", i, cell_coords(ndim*i)
        !write (6,*) "uval_current:", i, uval_current(i)
        label_current(i) = i
      enddo

      j = can_extent(0) !j表示与第0个单元相邻的结点索引在数组can_index中的起始位置.
      k = can_index(j)  !k表示一个结点的索引. 该结点在所有与第0个单元相邻的结点中排第0个.
      d_dx = 2.0d0 * dabs(cell_coords(0*ndim) - node_coords(k*ndim)); 

      return
      end
