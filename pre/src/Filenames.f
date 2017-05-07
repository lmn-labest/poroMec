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
c *                                                                    *
c **********************************************************************
      implicit none      
      include 'parallel.fi'
      character*80 NomeArqDados,NomeArqGerado
      character*20 StrExtensao
      integer      iPonto,TamanhoNome,NumArq,code
c ......................................................................
      if    (code .eq. 0) then
         write( StrExtensao, '( A )' ) '.geo'
      elseif(code .eq. 1) then                                
         write( StrExtensao, '( A,I4.4, A )' ) '_temp',NumArq,'.scl' 
      elseif(code .eq. 2) then
         write( StrExtensao, '( A )' ) '.res'        
      elseif(code .eq. 3) then
         write( StrExtensao, '( A,I4.4, A )' )'_temp.',NumArq, '.vtk' 
      elseif(code .eq. 4) then
         write( StrExtensao, '( A,I4.4, A )' ) '_flux',NumArq,'.vec' 
      elseif(code .eq. 5) then
         write( StrExtensao, '( A )' ) '.msh'        
      elseif(code .eq. 6) then
         write( StrExtensao, '( A,I4.4, A )' )'_flux.',NumArq, '.vtk' 
      elseif(code .eq. 7) then                                
         write( StrExtensao, '( A,I4.4, A )' ) '_disp',NumArq,'.scl'
      elseif(code .eq. 8) then
         write( StrExtensao, '( A )' ) '.vtk'
      elseif(code .eq. 9) then
         write( StrExtensao, '( A,I4.4, A )' )'_disp.',NumArq, '.vtk'
      elseif(code .eq. 10) then
         write( StrExtensao, '( A,I4.4, A )' ) '_stress',NumArq,'.scl'
      elseif(code .eq. 11) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        write( StrExtensao, '( 3A )' ) '_',trim(StrExtensao),'.dat'
       elseif(code .eq. 12) then
         write( StrExtensao, '( A,I4.4, A )' )'_stress.',NumArq, '.vtk'
      elseif(code .eq. 13) then
         write( StrExtensao, '( A,I4.4, A )' ) '_wave',NumArq,'.scl' 
      elseif(code .eq. 14) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
c        write( StrExtensao, '( 3A )' ) '_t_',trim(StrExtensao),'.txt'
        StrExtensao = '_t_'//trim(StrExtensao)//'.txt'
      elseif(code .eq. 15) then
         write( StrExtensao, '( A,I4.4, A )' )'_wave.',NumArq, '.vtk'
      elseif(code .eq. 17) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        write( StrExtensao, '( 3A )' ) '_log_',trim(StrExtensao),'.txt'
c ................................................................
c
c ... Codigos usados pelo PRE-processador:
c
      elseif(code .eq. 101) then                                
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao = '_'//trim(StrExtensao)//'.dat'
      elseif(code .eq. 102) then
        write(StrExtensao,'( I4 )') NumArq
        write(StrExtensao,'(A7,A)') '.epart.',trim(adjustl(StrExtensao))
      elseif(code .eq. 103) then
        write(StrExtensao,'( I4 )') NumArq
        write(StrExtensao,'(A7,A)') '.npart.',trim(adjustl(StrExtensao))
      elseif(code .eq. 104) then
        write(StrExtensao,'( I6 )') nprcs
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        if(ovlp)then
          write( StrExtensao, '( 3A )' )
     .      '_o',trim(StrExtensao),'_part.msh'
        else
          write( StrExtensao, '( 3A )' )
     .      '_n',trim(StrExtensao),'_part.msh'
        endif
      elseif(code .eq. 105) then
        write(StrExtensao,'( I6 )') nprcs
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        if(ovlp)then
          write( StrExtensao, '( 3A )' )
     .      '_o',trim(StrExtensao),'_part.res'
        else
          write( StrExtensao, '( 3A )' )
     .      '_n',trim(StrExtensao),'_part.res'
        endif
      elseif(code .eq. 106) then
        write(StrExtensao,'( I6 )') nprcs
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        if(ovlp)then
          StrExtensao = '_o_'//trim(StrExtensao)//'_part.vtk'
        else
          StrExtensao = '_n_'//trim(StrExtensao)//'_part.vtk'
        endif
      elseif(code .eq. 107) then
        write(StrExtensao,'( I6 )') NumArq 
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        if(ovlp)then
          StrExtensao = '_o_'//trim(StrExtensao)//'_my_part.vtk'
        else
          StrExtensao = '_n_'//trim(StrExtensao)//'_my_part.vtk'
        endif
      elseif(code .eq. 108) then                                
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao = '_map_'//trim(StrExtensao)//'.dat'
      elseif(code .eq. 109) then  
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao = '_pre_t_'//trim(StrExtensao)//'.txt'
      elseif(code .eq. 110) then
        write(StrExtensao,'( I6 )') nprcs
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
          StrExtensao = '_geo_'//trim(StrExtensao)//'_pre.vtk'
      elseif(code .eq. 111) then
        write(StrExtensao,'( I6 )') NumArq 
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        if(ovlp)then
          StrExtensao = '_o_'//trim(StrExtensao)//'_my_part.vtu'
        else
          StrExtensao = '_n_'//trim(StrExtensao)//'_my_part.vtu'
        endif
      elseif(code .eq. 116) then
        write(StrExtensao,'( I6 )') nprcs
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        if(ovlp)then
          StrExtensao = '_o_'//trim(StrExtensao)//'_part.vtu'
        else
          StrExtensao = '_n_'//trim(StrExtensao)//'_part.vtu'
        endif
      elseif(code .eq. 117) then
        write(StrExtensao,'( I6 )') nprcs
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
          StrExtensao = '_geo_'//trim(StrExtensao)//'_pre.vtu'
      elseif(code .eq. 118) then
        write(StrExtensao,'( I6 )') nprcs
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
          StrExtensao = '_geo_face_'//trim(StrExtensao)//'_pre.vtk'
      elseif(code .eq. 119) then
        write(StrExtensao,'( I6 )') nprcs
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
          StrExtensao = '_geo_face_'//trim(StrExtensao)//'_pre.vtu'
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
