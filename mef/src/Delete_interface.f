      subroutine del_interface(ix,numel_sint,numel,numel_i,nen,mat_i)
c **********************************************************************
c *                                                                    *
c *   DEL_INTERFACE:deleta os elementos de interface                   *
c *   -----                                                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ----------------------                                           *
c *                                                                    *
c *    ix(nen+1,numel) - conetividades nodais                          *
c *    x(ndm,nnode)    - coordenadas nodais                            *
c *    u(ndf,nnode)    - campo escalar ou vetorial                     *
c *    nnode  - numero de nos                                          *
c *    numel  - numero de elementos                                    *
c *    nquad4 - numero de quadrilateros                                *
c *    nprism6- numero de prismas                                      *
c *    nhexa8  - numero de hexaedros                                   *
c *    ndm - dimensao                                                  *
c *    nen - numero maximo de nos por elementos                        *
c *    ndf - dimensao de u                                             *
c *    nplot - arquivo de saida                                        *
c *    cont1 - controle de impressao vtk                               *
c *    code - codigo de instrucao:                                     *
c *       ---  Geometria ( valores < 0 )  ---                          *
c *         = -1 - geometria View3D                                    *
c *         = -2 - geometria Gid                                       *
c *         = -3 - geometria ParaView                                  *
c *       ---  Solucao ( valores > 0 )  ---                            *
c *         =  1 - Temperatura View3D                                  *
c *         =  2 - Temperatura Gid                                     *
c *         =  3 - Temperatura ParaView                                *
c *         =  4 - Fluxo de Calor View3D                               *
c *         =  5 - Fluxo de Calor Gid                                  *
c *         =  6 - Fluxo de Calor ParaView                             *
c *         =  7 - Deslocamentos View3D                                *
c *         =  8 - Deslocamentos Gid                                   *
c *         =  9 - Deslocamentos ParaView                              *
c *         =  10 - Tensoes View3D                                     *
c *         =  11 - Tensoes Gid                                        *
c *         =  12 - Tensoes ParaView                                   *
c *                                                                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   --------------------                                             *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'elementos.fi'
      integer ix(nen+1,*),numel,nen,mat_i
c      integer ix(21,*),nnode,numel,ndm,nen,ndf,nplot,code,istep
      integer numel_i,ndf
      integer i,j,k,n_inter,numel_sint(nen+1,numel)
c ======================================================================
      n_inter = 0
      do i = 1, numel
         if (ix(nen+1,i) .eq. 4)   then
            n_inter = n_inter + 1
         else 
            do j = 1, nen+1 
               numel_sint(j,i-n_inter) =ix(j,i)
            enddo
         endif
      enddo
      numel_i = numel - n_inter
c      print *, numel_i
      return
      end

          