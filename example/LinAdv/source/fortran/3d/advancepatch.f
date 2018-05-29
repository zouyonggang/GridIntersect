c**********************************************************************
c 计算通量.
c**********************************************************************
      subroutine advancepatch3d(velocity, 
     &                          dt, 
     &                          num_faces,
     &                          num_extend_cells, 
     &                          face_coords, 
     &                          cell_extend_coords, 
     &                          fac_extent, 
     &                          fac_index, 
     &                          uval,
     &                          flux)
c**********************************************************************
      implicit none
c**********************************************************************
c input parameters:
      double precision velocity, dt
      
      ! 网格片内部面的数目.
      integer num_faces
      ! 网格片内部及第1层影像区中的单元数之和.
      integer num_extend_cells

      ! <面, 单元>拓扑.
      integer fac_extent(0:num_faces)
      integer fac_index(0:fac_extent(num_faces)-1)

      ! 面中心坐标.
      double precision face_coords(0:num_faces*3-1)

      ! 单元中心坐标.
      double precision cell_extend_coords(0:num_extend_cells*3-1)

      ! 守恒量：中心量数据片
      double precision uval(0:num_extend_cells-1)

c output parameters:
      ! 通量：面心量数据片.
      double precision flux(0:num_faces-1)

c**********************************************************************
c     local parameters:
      integer i, idx1, idx2
      double precision tmp_coord
      double precision epsilon
      parameter (epsilon = 0.00000001)
c**********************************************************************

      do i = 0, num_faces-1
        flux(i) = 0.0d0
      enddo

      do i = 0, num_faces-1
        !用于识别的X, Y方向的边的辅助量.
        tmp_coord = face_coords(3*i) - epsilon

        if (fac_extent(i+1) - fac_extent(i) .eq. 2) then
          ! 如果与一个面相邻的单元有两个, 则说明这个面是内部面, 否则就是物理边界上的面.
          ! 取出与第i个面相邻的两个单元的索引.
          idx1 = fac_index(fac_extent(i))
          idx2 = fac_index(fac_extent(i)+1)
           
          ! 对于垂直于X坐标轴的面, 用其左边单元的中心量计算通量,
          ! 其它面上的通量值均为0.
          if (tmp_coord .gt. cell_extend_coords(3*idx1)) then 
            flux(i) = dt*uval(idx1)*velocity
          else if (tmp_coord .gt. cell_extend_coords(3*idx2)) then   
            flux(i) = dt*uval(idx2)*velocity
          endif

        ! 如果与一条边相邻的单元只有一个，则说明该边是物理边界上的边.
        else if (fac_extent(i+1) - fac_extent(i) .eq. 1) then
          idx1 = fac_index(fac_extent(i))
          if (tmp_coord .gt. cell_extend_coords(3*idx1)) then 
          ! X方向右边界上的面.
            flux(i) = dt * velocity * uval(idx1)
          endif
        else
          write (6,*) "Error:wrong number of cells neighbouring to face"
          stop
        endif
      enddo

      return
      end
               
