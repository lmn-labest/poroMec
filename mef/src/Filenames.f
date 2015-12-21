c*****************************Svn***************************************      
c*$Date: 2011-10-28 14:58:18 -0200 (Fri, 28 Oct 2011) $                 
c*$Rev: 956 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************      
      character*80 function name(NomeArqDados,NumArq,code)
c **********************************************************************
c *                                                                    *
c *   NAME: nomes de aquivos                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    NumArq       - numero do arquivo                                *
c *    code         - codigo de instrucao                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    NomeArqDados - nome do arquivo                                  *
c * 
c *    code:                                                           *
c *    0 -> geometria vtk                                              *
c *    1 -> geometrica vtu                                             *
c *    2 -> resultado vtk                                              *
c *    3 -> resultado vtu                                              *
c *    4 ->                                                            *
c *    5 ->                                                            *
c *    6 ->                                                            *
c *    7 ->                                                            *
c *    8 ->                                                            *
c *    9 ->                                                            *
c *   10 ->                                                            *
c *   12 ->                                                            *
c *   13 ->                                                            *
c *   14 -> arquivo de tempos                                          *
c *   15 -> arquivo do log do solver                                   *
c *   30 -> arquivos de resultados por nos                             *
c *   31 -> arquivos de tensoes total por nos                          *
c *   32 -> arquivos de tensoes efetivas de terzaghi por nos           *
c *   33 -> arquivos de tensoes efetivas de biot por nos               *
c *   34 -> arquivos de fluxo de darcy por nos                         *
c **********************************************************************
      implicit none      
      include 'parallel.fi'
      character*80 NomeArqDados,NomeArqGerado
      character*30 StrExtensao
      integer      iPonto,TamanhoNome,NumArq,code
c ......................................................................
c
c ... geom .vtk
      if    (code .eq. 0) then
        StrExtensao='_geo.vtk'
c ... geom .vtu
      else if    (code .eq. 1) then
        StrExtensao='_geo.vtu'
c ... res .vtk
      else if    (code .eq. 2) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)  
        StrExtensao='_istep_'//trim(StrExtensao)//'.vtk'
c ... res .vtu
      else if    (code .eq. 3) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)  
        StrExtensao='_istep_'//trim(StrExtensao)//'.vtu'
c ... log de tempos
      elseif(code .eq. 14) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_t_'//trim(StrExtensao)//'.txt'
c ... log solver
      elseif(code .eq. 15) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_log_'//trim(StrExtensao)//'.txt'
c ... pu node file 
      elseif(code .eq. 30) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_up_node_'//trim(StrExtensao)//'.txt'
c ... stress node file 
      elseif(code .eq. 31) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_stress_node_'//trim(StrExtensao)//'.txt'
c ... terzaghi stress node file 
      elseif(code .eq. 32) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_stressE_node_'//trim(StrExtensao)//'.txt'
c ... biot stress node file 
      elseif(code .eq. 33) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_stressB_node_'//trim(StrExtensao)//'.txt'
c ... flux node file 
      elseif(code .eq. 34) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_flux_node_'//trim(StrExtensao)//'.txt'
c ... mtx 
      elseif(code .eq. 50) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_'//trim(StrExtensao)
c ... kuu mtx 
      elseif(code .eq. 51) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_kuu_'//trim(StrExtensao)
c ... kpp mtx 
      elseif(code .eq. 52) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_kpp_'//trim(StrExtensao)
      endif
c ......................................................................      
      TamanhoNome = INDEX( NomeArqDados, ' '  )
      iPonto = INDEX( NomeArqdados, '.' )
      if( iPonto .EQ. 0 ) then
          NomeArqGerado = NomeArqDados(1:TamanhoNome-1) // StrExtensao
      else
          NomeArqGerado = NomeArqDados(1:iPonto-1) // StrExtensao
      endif
      name = NomeArqGerado
c ......................................................................      
      return
      end
