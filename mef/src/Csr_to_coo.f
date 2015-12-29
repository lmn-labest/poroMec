c********************************************************************* 
c* CSR_TO_COO_PM : conveter do formato CSR para COO                  * 
c*-------------------------------------------------------------------* 
c* Parametros de entrada:                                            * 
c*-------------------------------------------------------------------* 
c* lin   -> indefinido                                               * 
c* col   -> indefinido                                               * 
c* val   -> indefinido                                               * 
c* ia    -> vetor CSR                                                * 
c* ja    -> vetor CSR                                                * 
c* ad    -> matrix de coeficientes                                   * 
c* al    -> matrix de coeficientes                                   * 
c* neq   -> numero de equacoes                                       * 
c* nad   -> numero de termos nao nulos                               * 
c* bin   -> matriz binaria                                           * 
c*-------------------------------------------------------------------* 
c* Parametros de saida:                                              * 
c*-------------------------------------------------------------------* 
c* lin -> numero da linha                                            * 
c* col -> numero da coluna                                           * 
c* val -> valor                                                      * 
c*-------------------------------------------------------------------* 
c* OBS:                                                              * 
c*-------------------------------------------------------------------* 
c********************************************************************* 
      subroutine csr_to_coo_pm(lin   ,col,val
     .                        ,ia    ,ja
     .                        ,al    ,ad 
     .                        ,neq   ,nad,bin)
      implicit none
      integer lin(*),col(*)
      real*8  val(*)
      integer ia(*),ja(*)
      real*8 al(*),ad(*)   
      integer neq ,nad
      logical unsym,bin,pre_diag
      integer nl,nc,kk
c     pre_diag = .true.
      pre_diag = .false.
c ... CSRC poro mecanico (ad-diagonal principal;al-parte inferior)
      kk = 0
      do nl = 1, neq
        kk     = kk + 1
        lin(kk) = nl
        col(kk) = nl
        if(bin) then
          val(kk) = 1.0
        else
          if(pre_diag) then
            val(kk) = 1.d0
          else
            val(kk) = ad(nl)
          endif
        endif
        do nc = ia(nl), ia(nl+1) - 1
          kk = kk + 1
          lin(kk)  = nl
          col(kk)  = ja(nc)
          if(bin) then
            val(kk) = 1.0
          else
            if(pre_diag) then
               val(kk) = al(nc)/ad(nl)
            else
              val(kk) = al(nc)
            endif
          endif  
          kk       = kk + 1
          lin(kk)  = ja(nc)
          col(kk)  = nl        
          if(bin) then
            val(kk) = 1.0
          else
            if(pre_diag) then
              val(kk) = al(nc)/ad(ja(nc))
            else
              val(kk) = al(nc)
            endif  
          endif
        enddo
      enddo
c .....................................................................
c
c ...
      return
      end
c *********************************************************************
c
c********************************************************************* 
c* CSR_BOCK_TO_COO : conveter do formato CSR(Kuu, kpp e kup)         * 
c * para COO(Kuu,Kpp, e K)                                           * 
c*-------------------------------------------------------------------* 
c* Parametros de entrada:                                            * 
c*-------------------------------------------------------------------* 
c* lin   -> indefinido                                               * 
c* col   -> indefinido                                               * 
c* val   -> indefinido                                               * 
c* ia    -> vetor CSR                                                * 
c* ja    -> vetor CSR                                                * 
c* au    -> matrix de coeficientes                                   * 
c* ad    -> matrix de coeficientes                                   * 
c* al    -> matrix de coeficientes                                   * 
c* neq   -> numero de equacoes                                       * 
c* nad   -> numero de termos nao nulos                               * 
c* code  -> 1 - CSRD; 2 - CSR; 3 - CSRC                              * 
c* unsym -> matriz nao simetrica                                     * 
c* bin   -> matriz binaria                                           * 
c*-------------------------------------------------------------------* 
c* Parametros de saida:                                              * 
c*-------------------------------------------------------------------* 
c* lin -> numero da linha                                            * 
c* col -> numero da coluna                                           * 
c* val -> valor                                                      * 
c*-------------------------------------------------------------------* 
c* OBS:                                                              * 
c*-------------------------------------------------------------------* 
c********************************************************************* 
      subroutine csr_block_to_coo(lin   ,col  ,val
     .                           ,linuu ,coluu,valuu
     .                           ,linpp ,colpp,valpp
     .                           ,ia    ,ja
     .                           ,al    ,ad 
     .                           ,neq   ,nequ ,nad,nadpu
     .                           ,bin)
      implicit none
      integer lin(*),col(*),linuu(*),coluu(*),linpp(*),colpp(*)
      real*8  val(*),valuu(*),valpp(*)
      integer ia(*),ja(*)
      real*8 al(*),ad(*)   
      integer neq,nequ,nad,nadpu
      logical bin,pre_diag
      integer nl,nc,kk,ii
      integer desloc_ia,desloc_ja
      pre_diag = .false.
c ... matriz completa
      kk = 0
c ... Kuu e Kpp
      do nl = 1, neq
        kk     = kk + 1
        lin(kk) = nl
        col(kk) = nl
        if(bin) then
          val(kk) = 1.0
        else
          if(pre_diag) then
            val(kk) = 1.d0
          else
            val(kk) = ad(nl)
          endif
        endif
        do nc = ia(nl), ia(nl+1) - 1
          kk = kk + 1
          lin(kk)  = nl
          col(kk)  = ja(nc)
          if(bin) then
            val(kk) = 1.0
          else
            if(pre_diag) then
              val(kk) = al(nc)/ad(nl)
            else
              val(kk) = al(nc)
            endif  
          endif
          kk = kk + 1
          lin(kk)  = ja(nc)
          col(kk)  = nl     
          if(bin) then
            val(kk) = 1.0
          else
            if(pre_diag) then
              val(kk) = al(nc)/ad(ja(nc))
            else
              val(kk) = al(nc)
            endif  
          endif
        enddo
      enddo
c .....................................................................
c
c ... kpu
      desloc_ia = neq+1
      desloc_ja = nad
      do nl = 1, neq - nequ
        do nc = ia(desloc_ia+nl), ia(desloc_ia+nl+1) - 1
          kk       = kk + 1
          lin(kk)  = nl + nequ
          col(kk)  = ja(desloc_ja+nc)
          if(bin) then
            val(kk) = 1.0
          else
            if(pre_diag) then
              val(kk) = al(desloc_ja+nc)/ad(nl+nequ)
            else
              val(kk) = al(desloc_ja+nc)
            endif  
          endif
c
          kk       = kk + 1
          lin(kk)  = ja(desloc_ja+nc)
          col(kk)  = nl + nequ
          if(bin) then
            val(kk) = 1.0
          else
            if(pre_diag) then
              val(kk) = -al(desloc_ja+nc)/ad(ja(desloc_ja+nc))
            else
              val(kk) = -al(desloc_ja+nc)
            endif  
          endif
        enddo
      enddo 
c .....................................................................
c
c ... matriz Kuu      
      kk = 0
c ... Kuu      
      do nl = 1, nequ
        kk        = kk + 1
        linuu(kk) = nl
        coluu(kk) = nl
        if(bin) then
          valuu(kk) = 1.0
        else
          if(pre_diag) then
            valuu(kk) = 1.d0
          else
            valuu(kk) = ad(nl)
          endif
        endif
        do nc = ia(nl), ia(nl+1) - 1
          kk         = kk + 1
          linuu(kk)  = nl
          coluu(kk)  = ja(nc)
          if(bin) then
            valuu(kk) = 1.0
          else
            if(pre_diag) then
              valuu(kk) = al(nc)/ad(nl)
            else
              valuu(kk) = al(nc)
            endif  
          endif
          kk = kk + 1
          linuu(kk)  = ja(nc)
          coluu(kk)  = nl     
          if(bin) then
            valuu(kk) = 1.0
          else
            if(pre_diag) then
              valuu(kk) = al(nc)/ad(ja(nc))
            else
              valuu(kk) = al(nc)
            endif  
          endif
        enddo
      enddo
c .....................................................................
c
c ... matriz Kpp      
      kk = 0
c ... Kuu      
      do nl = 1, neq-nequ
        ii        = nequ + nl
        kk        = kk + 1
        linpp(kk) = nl
        colpp(kk) = nl
        if(bin) then
          valpp(kk) = 1.0
        else
          if(pre_diag) then
            valpp(kk) = 1.d0
          else
            valpp(kk) = ad(ii)
          endif
        endif
        do nc = ia(ii), ia(ii+1) - 1
          kk         = kk + 1
          linpp(kk)  = nl
          colpp(kk)  = ja(nc) - nequ
          if(bin) then
            valpp(kk) = 1.0
          else
            if(pre_diag) then
              valpp(kk) = al(nc)/ad(ii)
            else
              valpp(kk) = al(nc)
            endif  
          endif
          kk = kk + 1
          linpp(kk)  = ja(nc) - nequ
          colpp(kk)  = nl     
          if(bin) then
            valpp(kk) = 1.0
          else
            if(pre_diag) then
              valpp(kk) = al(nc)/ad(ja(nc)-nequ)
            else
              valpp(kk) = al(nc)
            endif  
          endif
        enddo
      enddo
c .....................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************
c* WRITECOO : escrita do sistem linear no formato COO                * 
c*-------------------------------------------------------------------* 
c* Parametros de entrada:                                            * 
c*-------------------------------------------------------------------* 
c* lin     -> numero da linha COO                                    * 
c* col     -> numero da coluna COO                                   * 
c* val     -> valor COO                                              * 
c* f       -> vetor de forcas                                        * 
c* neq     -> numero de equacoes                                     * 
c* nnz     -> numero totoal de elementos nao nulos                   * 
c* nad     -> numero de termos nao nulos                             * 
c* nameout ->                                                        * 
c* fileout -> numero do arquivo de saida                             * 
c*  flag   -> escreve o vetor de forcas                              * 
c* unsym   -> matriz nao simetrica                                   * 
c*-------------------------------------------------------------------* 
c* Parametros de saida:                                              * 
c*-------------------------------------------------------------------* 
c*-------------------------------------------------------------------* 
c* OBS:                                                              * 
c*-------------------------------------------------------------------* 
c********************************************************************* 
      subroutine write_coo(lin    ,col    ,val
     .                    ,f      ,neq    ,nnz
     .                    ,prefixo,fileOut
     .                    ,flag   ,unsym  )
      implicit none
      integer lin(*),col(*),nnz,i,neq
      real*8  val(*),f(*)
      integer fileOut
      logical flag,unsym
      character*80  prefixo
      character*100 nameOut
      character*1024 str
c
      nameOut = trim(prefixo) //'.mtx'
c ... matriz de coeficientes
      open(fileOut , action='write', file=nameOut)
      if(unsym) then
        str = '%%MatrixMarket matrix coordinate real general'
      else
        str = '%%MatrixMarket matrix coordinate real symmetric'
      endif                                                      
      write(fileOut,'(a)')trim(str)
      write(fileOut,'(3i9)')neq,neq,nnz
      do i = 1, nnz
        write(fileOut,'(2i9,e24.16)')lin(i),col(i),val(i)
      enddo
      close(fileOut)
c .....................................................................
c
c ...
c .....................................................................
c
c ... vetor de forca
      if( flag ) then
        nameOut = trim(prefixo) //'_b.mtx'
        open(fileOut , action='write', file=nameOut)
        str ='%%MatrixMarket matrix array real general'
        write(fileOut,'(a)')trim(str)
        write(fileOut,'(2i9)')neq,1
        do i = 1, neq
          write(fileOut,'(e24.16)') f(i)
        enddo
        close(fileOut)
      endif
c .....................................................................
      return
      end
c *********************************************************************
