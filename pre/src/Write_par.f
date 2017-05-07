c *********************************************************************
c * MYTEMP    : estatisca de tempo de excucao                         *
c * ----------------------------------------------------------------- *
c * Parametro de entrada :                                            *
c * ----------------------------------------------------------------- *
c * trmetis       - tempo na leitura da malha do metis                *
c * tcomunica     - tempo nas comunicacao de Mpi                      *
c * tpart         - tempo na geracao dos mapas                        *
c * tmetis        - tempo no metis                                    *
c * trmesh        - tempo no leitura para o mefpar                    *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * p (el)           - matererias por elementos locais                *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine mytemp(MPIW,filein,nnos,rank,nprcs,naux)
c ===
      implicit none
      include 'time.fi'
      real*8 total
      integer MPIW
      character*80 filein,name
      character*100 fileout
      integer nnos,rank
      integer nprcs
      integer naux
      integer ierr
c =====================================================================
c
c ===
#ifdef MPI
      real*8 tc,tm,tp
c ... contribuicoes de cada processo
      tc        = tcomunica
      tp        = tpart
      tcomunica = 0.d0
      tmetis    = 0.d0
      tpart     = 0.d0
c .....................................................................
c
c ...
      call Mpi_temp(MPIW,tcomunica,tmetis,tpart,tc,tm,tp,ierr)
c .....................................................................      
#endif
c =====================================================================
c
c ===
      if( rank .eq. 0 ) then
        fileout = name(filein,nprcs,109)
        open(naux,file=fileout)
c =====================================================================
c
c ===
        total = trmetis + tcomunica/nnos + tmetis/nnos + tpart/nnos
     .        + trmesh 
        write(naux,'(A)')'Tempo'
        write(naux,'(A,f30.5)')'trmetis  ',trmetis
        write(naux,'(A,f30.5)')'tcomunica',tcomunica/nnos
        write(naux,'(A,f30.5)')'tpart    ',tpart/nnos
        write(naux,'(A,f30.5)')'trmesh   ',trmesh
        write(naux,'(A,f30.5)')'Total    ',total  
c =====================================================================
c
c ===
        ierr = 0
        close(naux)
      endif  
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * LOCALCOOR : gera as coordenadas locais                            *
c * ----------------------------------------------------------------- *
c * Parametro de entrada :                                            *
c * ----------------------------------------------------------------- *
c * x(ndm,nnode)     - coordenadas locais                             *
c * noLG(my_nnode)   - map local -> global de nos                     *
c * ndm              - numero de dimensoes                            *
c * my_node          - numero de nos locais                           *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * xl(ndm,my_nnode) - coordenadas locais                             *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine localcoor(x,xl,noLG,my_nnode,ndm)
c ===
      implicit none
c ... malha
      real*8  x(ndm,*)
      integer ndm
c ... malha local
      real*8  xl(ndm,*)
      integer my_nnode
c ... particionamento
      integer noLG(*)
c ... auxiliar
      integer i,k
c =====================================================================
c
c ===
      do i = 1, my_nnode
        k = noLG(i)
        call copy_coor(x,xl,i,k,ndm)
      enddo
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * LOCALEL : gera as conectividades locais                           *
c * ----------------------------------------------------------------- *
c * Parametro de entrada :                                            *
c * ----------------------------------------------------------------- *
c * el(nen,numel)  - conectividodes global                            *
c * noGL(nnode)    - map global -> local de nos                       *
c * elLG(my_nnode) - map local -> global dos elementos                *
c * ep(numel)      - particionamento por elemento                     *
c * my_numel       - numero de elementos locais                       *
c * nen            - numero de nos por elementos                      *
c * iws            - escreve o particionamento(false ou true)         *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * lel(nen+1,my_numel) - elementos locais                            *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine localel(el,lel,elLG,noGL,my_numel,nen)
c ===
      implicit none
c ... malha
      integer el(nen+1,*)
      integer nen
c ... malha local
      integer lel(nen+1,*)
      integer my_numel
c ... particionamento
      integer elLG(*),noGL(*)
c ... auxiliar
      integer i,k
c =====================================================================
c
c ===
      do i = 1, my_numel
        k = elLG(i)
        call copy_elm(el,lel,i,k,nen)
        call elmGL(noGL,lel,i,nen)
      enddo
c ====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * LOCAL_NO_PMEC : gera as propriedades locais do poro mec           *
c * ----------------------------------------------------------------- *
c * Parametro de entrada :                                            *
c * ----------------------------------------------------------------- *
c * el(nen,numel)  - conectividodes global                            *
c * u(nnode)       - proprideade global                               * 
c * lu             - nao definido                                     * 
c * noGL(nnode)    - map global -> local de nos                       *
c * elLG(my_nnode) - map local -> global dos elementos                *
c * ep(numel)      - particionamento por elemento                     *
c * my_numel       - numero de elementos locais                       *
c * maxnov         - numero de nos de vertices                        *
c * maxno          - numero de nos por elementos                      *
c * ndf            - graus de liberade em u                           *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * lu (my_numel)  - propriedade nodal locais                         *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine local_no_pmec(el      ,u     ,lu   ,elLG,noGL
     .                        ,my_numel,maxnov,maxno,ndf)
c ===
      implicit none
c ... malha
      real*8 u(ndf,*)
      integer el(maxno+1,*)
      integer maxno,maxnov,ndf
c ... malha local
      real*8 lu(*)
      integer my_numel
c ... particionamento
      integer elLG(*),noGL(*)
c ... auxiliar
      integer i,j,k,no
c =====================================================================
c
c ===
      do i = 1, my_numel
        k = elLG(i)
        do j = 1, maxnov
          no = el(j,k)  
          if( no .ne. 0 ) then
            lu(noGL(no)) = u(ndf,no)
          endif
        enddo 
      enddo
c ====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * ELMGL   : gera as conectividades locais                           *
c * ----------------------------------------------------------------- *
c * Parametro de entrada :                                            *
c * ----------------------------------------------------------------- *
c * noGL(nnode)    - map global -> local de nos                       *
c * lel(nen+1,*)   - conectividas globais dos elementos locais        *
c * nlel           - numero do elemento locais                        *
c * nen            - numero de nos por elementos                      *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * lel(nen+1,*)   - conectividas locais do elementos locais          *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine elmGL(noGL,lel,nlel,nen)
c ===
      implicit none
      integer noGL(*),lel(nen+1,*)
      integer nlel,nen
      integer no, lno
      integer j
c =====================================================================
c
c ===
      do j = 1 , nen
        no  = lel(j,nlel)
        lno = noGL(no)
        lel(j,nlel) = lno
      enddo
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * COPY_ELM: copia os conectividades globais para o arranjo local    *
c * ----------------------------------------------------------------- *
c * Parametro de entrada :                                            *
c * ----------------------------------------------------------------- *
c * el(nen,numel)       - conectividodes globais                      *
c * lel(nen+1,my_numel) - conectividades locais                       *
c * nlel                - numero do elemento local                    *
c * nel                 - numero do elemento global                   *
c * ep(numel)           - particionamento por elemento                *
c * nen                 - numero de nos por elementos                 *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * lel(nen+1,my_numel) - elementos locais                            *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine copy_elm(el,lel,nlel,nel,nen)
c ===
      implicit none
      integer el(nen+1,*),lel(nen+1,*)
      integer nel,nen,nlel
      integer j
c =====================================================================
c
c ===
      do j = 1, nen + 1
        lel(j,nlel) = el(j,nel)
      enddo
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * COPY_COOR : copia os coordenadas  globais para o arranjo local    *
c * ----------------------------------------------------------------- *
c * Parametro de entrada :                                            *
c * ----------------------------------------------------------------- *
c * x(ndm,nnode)  - coordenadas globais                               *
c * xl(ndm,nnode) - coordenadas locais                                *
c * nlno          - numero do ponto local                             *
c * nno           - numero do ponto global                            *
c * ndm           - numero de dimensoes                               *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * xl(ndm,my_nnode) - pontos locais                                  *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine copy_coor(x,lx,nlno,nno,ndm)
c ===
      implicit none
      real*8  x(ndm,*),lx(ndm,*)
      integer nno,ndm,nlno
      integer j
c =====================================================================
c
c ===
       do j = 1, ndm
          lx(j,nlno) = x(j,nno)
       enddo
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
c *********************************************************************
c * COPY_MAT  : copia os elementos locais                             *
c * ----------------------------------------------------------------- *
c * Parametro de entrada :                                            *
c * ----------------------------------------------------------------- *
c * p(ndm,nnode)  - materias por lementos locais                      *
c * lel(*,*)      - conectividades locais(caracteristica do mef)      *
c * my_numel      - numero de elementos locais                        *
c * nen           - numero de nos no elemento                         *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * p (el)           - matererias por elementos locais                *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine copy_mat(p,lel,my_numel,nen)
c ===
      implicit none
      integer p(*)
      integer lel(nen+1,*)
      integer my_numel,nen
      integer i
c =====================================================================
c
c ===
      do i = 1, my_numel
        p(i) = lel(nen+1,i)
      enddo  
c =====================================================================
c
c ===
      return
      end
c =====================================================================
c *********************************************************************
c
