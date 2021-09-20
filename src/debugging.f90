module debugging
    implicit none
    !
    interface mat_dump
        module procedure mat_dump_mat,mat_dump_mat_2_file,mat_dump_vec, mat_dump_vec_2_file
    end interface

    CONTAINS

!statusbar(
!            now (int)='position up to now',
!            tot (int)='total length to achieve, e.g. the end value of the do loop',
!            barlen (int,opt.)='number of character one wants'
!          )
     subroutine statusbar(now,tot,barlen)
        integer*4, intent(in)                :: tot,now
        integer*4, intent(in), optional        :: barlen
        !character*1, parameter                :: symbol="#",void=""
        character*150                        :: bar
        integer*4                            :: ii, entries, barstatuslen=20
        
        if(present(barlen)) barstatuslen=barlen
        entries = barstatuslen * now / tot
        
        bar="["
        do ii = 1,barstatuslen
            if (ii<entries)then
                bar=trim(bar)// "="
            elseif (ii==entries)then
                bar=trim(bar)// ">"
            else
                bar=trim(bar)// "--"
            endif
            if (ii == barstatuslen) then
                bar=trim(bar)// "]"
            endif
        enddo
        
        write(*,'(I4,2A)')100 * now / tot,"%",bar
        return
    endsubroutine statusbar
    
!date_time(
!            )
!scrive sul terminale la data e l'ora del momento in cui avviene la chiamata
    subroutine date_time()    
        character(8)            :: date
        character(10)           :: time
        character(5)            :: zone
        integer,dimension(8)    :: values
        ! using keyword arguments
        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)
        call date_and_time(VALUES=values)
        !print '(a,2x,a,2x,a)', date, time, zone
        !write(*,*) '*****************'
        write(*, '("***  ",I0,"/",I0,"/",I0,"  ",I0,":",I0,":",I0,"  ***")') values(3),values(2),values(1),&
                                                                            values(5),values(6),values(7)
        !write(*,*) '*****************'
        return
    end subroutine

!checkpoint( 
!             DEBUG='variabile logica',
!             TEXTOK (opt) = 'testo se il debu è positivo',
!             TEXT='testo da essere stampato'
!             var(opt)='variabile opzionale da stampare',
!             vec(opt),mat (opt) = 'come prima, ma rispettivamente array, matrice'
!             flag(int, opt) = 'flag=1->program continue, otherwise stops (default=0)'
!            )
    subroutine checkpoint(DEBUG,TEXTOK,TEXT,var,vec,mat,flag)
        implicit none
        !inputs
        LOGICAL, intent(IN)                             :: DEBUG
        CHARACTER (*), INTENT(IN)                        :: TEXT
        CHARACTER (*), INTENT(IN),optional                :: TEXTOK
        class(*),INTENT(IN), optional                    :: var
        class(*),INTENT(IN), dimension(:), optional        :: vec
        class(*),INTENT(IN), dimension(:,:), optional    :: mat
        integer,intent(in),optional                        :: flag
        !        
        IF (DEBUG) THEN
            if(present(TEXTOK))then
                WRITE(*,*)'* [DEBUG OK]****************************************************'
                WRITE(*,*)'* [MESSAGE] "',TEXTOK,'"'
                WRITE(*,*) '***************************************************************'
            endif
        ELSE 
            WRITE(*,*)'* [DEBUG ERR]**************************************************'
            WRITE(*,*)'* [ERROR MSG] "',TEXT,'"'
            if(present(var))then
                select type (var)
                    type is (integer(2))
                        write(*,'(" * [ERROR VAR]  ",I0," (integer*2)")') var 
                    type is (integer(4))
                        write(*,'(" * [ERROR VAR]  ",I0," (integer*4)")') var
                    type is (integer(8))
                        write(*,'(" * [ERROR VAR]  ",I0," (integer*8)")') var
                    type is (real(4))
                        write(*,'(" * [ERROR VAR]  ",F20.15," (real*4)")') var
                    type is (real(8))
                        write(*,'(" * [ERROR VAR]  ",F20.15," (real*8)")') var
                    type is (real(16))
                        write(*,'(" * [ERROR VAR]  ",F20.15," (real*16)")') var
                    type is (complex(8))
                        write(*,'(" * [ERROR VAR]  ",F20.15," (complex*8)")') var
                    type is (complex(16))
                        write(*,'(" * [ERROR VAR]  ",F20.15," (complex*16)")') var
                    type is (logical)
                        write(*,'(" * [ERROR VAR]  ",L," (logical)")') var
                    type is (character(*))
                        write(*,'(" * [ERROR VAR]  ",A," (character)")') var
                endselect 
            endif
            if(present(vec))then
                call mat_dump(vec=vec,text='* [IN DEBUG]: the input array is')
            endif
            if(present(mat))then
                call mat_dump(mat=mat,text='* [IN DEBUG]: the input matrix is')
            endif
            if(present(flag))then
                select case (flag)
                    case(0)
                        WRITE(*,*) '* [PRG STATUS] --- PROGRAM HAS STOPPED ---  '
                        stop
                    case(1)
                        WRITE(*,*) '* [PRG STATUS] --- PROGRAM CONTINUE ---  '
                    case default
                        WRITE(*,*) '* [PRG STATUS] --- PROGRAM HAS STOPPED ---  '
                        stop
                endselect
            else
                WRITE(*,*) '* [PRG STATUS] --- PROGRAM HAS STOPPED --- '
                WRITE(*,*) '***************************************************************'
                stop
            endif
            WRITE(*,*) '***************************************************************'
        ENDIF
        return 
    endsubroutine checkpoint

!mat_dump_mat(    mat='input matrix'
!                 text='optional text to put with the name of the matrix'
!                 opt='print format of matrix'
!                    opt=1-> (-8.960E+00,  9.500E+00)
!                   opt=2-> (    -8.960,      9.500)
!                   opt=3-> ( 9.140,  0.000)
!                       opt=4-> ( 0.5946107375736628, 0.6855384522353541)
!                flag='integer: 1-> no info flag, 0 or nothing-> info flag'
!              )
    subroutine mat_dump_mat(mat,text,opt,flag)
        implicit none
        !inputs
        class(*), dimension(:,:), intent(in) :: mat
        character(*), intent(in), optional :: text
        integer, intent(in), optional :: opt
        integer, intent(in), optional :: flag
        !system
        integer*4 ii,jj,N,N2
        character(50) fmtStringCmplx,fmtString,fmtStringInt,typeof
        !
        if (present(opt))then!può essere aumentato in futuro
           select case (opt)
                case(1)
                    fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
                    fmtString = '(9999("(",ES12.4,")",:,1X))'
                case(2)
                    fmtStringCmplx = '(9999("(",F10.3,",",1X,F10.3,")",:,1X))'
                    fmtString = '(9999("(",F10.3,")",:,1X))'
                case(3)
                    fmtStringCmplx = '(9999("(",F3.1,",",1X,F3.1,")",:,1X))'                      
                    fmtString = '(9999("(",F3.1,")",:,1X))'
                case(4)
                    fmtStringCmplx = '(9999("(",F20.16,",",1X,F20.16,")",:,1X))'                      
                    fmtString = '(9999("(",F20.16,")",:,1X))'
                case(5)!mathematica
                    fmtStringCmplx = '("{",9999(F20.16,"+ I(",1X,F20.16,"),",:,1X),"}")'                      
                    fmtString = '(9999("{",F20.16,",",:,1X),"}")'
                case default
                    fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
                    fmtString = '(9999("(",ES12.4,")",:,1X))'
           endselect
        else
            fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
            fmtString = '(9999("(",ES12.4,")",:,1X))'
        endif
        fmtStringInt = '(9999("(",I0,")",:,1X))'                  

        !print*,''
        print*,text
        typeof='unknown'
        
        N=size(mat,1)
        N2=size(mat,2)
        select type(mat)
            type is (complex(16))
                typeof='complex*16'
                do ii=1,N
                   write(*,fmt=fmtStringCmplx) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (complex(8))
                typeof='complex*8'
                do ii=1,N
                   write(*,fmt=fmtStringCmplx) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (complex(4))
                typeof='complex*4'
                do ii=1,N
                   write(*,fmt=fmtStringCmplx) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (integer(2))
                typeof='integer*2'
                do ii=1,N
                   write(*,fmt=fmtStringInt) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (integer(4))
                typeof='integer*4'
                do ii=1,N
                   write(*,fmt=fmtStringInt) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (real(4))
                typeof='real*4'
                do ii=1,N
                    write(*,fmtString) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (real(8))
                typeof='real*8'
                do ii=1,N
                    write(*,fmtString) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (real(16))
                typeof='real*16'
                do ii=1,N
                   write(*,fmt=fmtString) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (character(len=*))
                typeof='character'
                do ii=1,N
                    write(*,fmt='(9999("(",A,")",:,1X))') (mat(ii,jj),jj=1,size(mat,2))
                enddo
            type is (logical)
                typeof='logical'
                do ii=1,N
                    write(*,fmt='(9999("(",L,")",:,1X))') (mat(ii,jj),jj=1,size(mat,2))
                enddo
            end select 
        if(present(flag))then
            select case (flag)
            !no info message
                case(1)                   
                    print*,''
                case default
                    print*,'attenzione: flag/=1'
            endselect
        else
            write(*,'(" [matrix details: ",A,", dim: ",I0,"x",I0,"]")') trim(typeof),N,N2
        endif
        return
    endsubroutine mat_dump_mat

    subroutine mat_dump_mat_2_file(mat,unit_file,text,opt,flag)
        implicit none
        !inputs
        class(*), dimension(:,:), intent(in) :: mat
        integer, intent(in)           :: unit_file
        character(*), intent(in), optional :: text
        integer, intent(in), optional :: opt
        integer, intent(in), optional :: flag
        !system
        integer*4 ii,jj,N
        character(50) fmtStringCmplx,fmtString,fmtStringInt,typeof
        !
        if (present(opt))then!può essere aumentato in futuro
           select case (opt)
                case(1)
                    fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
                    fmtString = '(9999("(",ES12.4,")",:,1X))'
                case(2)
                    fmtStringCmplx = '(9999("(",F10.3,",",1X,F10.3,")",:,1X))'
                    fmtString = '(9999("(",F10.3,")",:,1X))'
                case(3)
                    fmtStringCmplx = '(9999("(",F3.1,",",1X,F3.1,")",:,1X))'                      
                    fmtString = '(9999("(",F3.1,")",:,1X))'
                case(4)
                    fmtStringCmplx = '(9999("(",F20.16,",",1X,F20.16,")",:,1X))'                      
                    fmtString = '(9999("(",F20.16,")",:,1X))'
                case(5)!mathematica
                    fmtStringCmplx = '("{",9999(F20.16,"+ I(",1X,F20.16,"),",:,1X),"}")'                      
                    fmtString = '(9999("{",F20.16,",",:,1X),"}")'
                case default
                    fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
                    fmtString = '(9999("(",ES12.4,")",:,1X))'
           endselect
        else
            fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
            fmtString = '(9999("(",ES12.4,")",:,1X))'
        endif
        fmtStringInt = '(9999("(",I0,")",:,1X))'                  

        !print*,''
        write(unit_file,*)text
        typeof='unknown'
        
        N=size(mat,1)
        select type(mat)
            type is (complex(16))
                typeof='complex*16'
                do ii=1,N
                   write(unit_file,fmt=fmtStringCmplx) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (complex(8))
                typeof='complex*8'
                do ii=1,N
                   write(unit_file,fmt=fmtStringCmplx) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (complex(4))
                typeof='complex*4'
                do ii=1,N
                   write(unit_file,fmt=fmtStringCmplx) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (integer(2))
                typeof='integer*2'
                do ii=1,N
                   write(unit_file,fmt=fmtStringInt) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (integer(4))
                typeof='integer*4'
                do ii=1,N
                   write(unit_file,fmt=fmtStringInt) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (real(4))
                typeof='real*4'
                do ii=1,N
                    write(unit_file,fmtString) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (real(8))
                typeof='real*8'
                do ii=1,N
                    write(unit_file,fmtString) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (real(16))
                typeof='real*16'
                do ii=1,N
                   write(unit_file,fmt=fmtString) ( mat( ii, jj ), jj = 1, size(mat,2) )
                enddo
            type is (character(len=*))
                typeof='character'
                do ii=1,N
                    write(unit_file,fmt='(9999("(",A,")",:,1X))') (mat(ii,jj),jj=1,size(mat,2))
                enddo
            type is (logical)
                typeof='logical'
                do ii=1,N
                    write(unit_file,fmt='(9999("(",L,")",:,1X))') (mat(ii,jj),jj=1,size(mat,2))
                enddo
            end select 
        if(present(flag))then
            select case (flag)
            !no info message
                case(1)                   
                    print*,''
                case default
                    print*,'attenzione: flag/=1'
            endselect
        else
            write(unit_file,'(" [matrix details: ",A,", dim: ",I0,"x",I0,"]")') trim(typeof),N,size(mat,2)
        endif
        return
    endsubroutine mat_dump_mat_2_file
!mat_dump_vec(
!                 vec='input vec'
!                 text='optional text to put with the name of the matrix'
!                 opt='print format of matrix'
!                     opt=1-> (-8.960E+00,  9.500E+00)
!                      opt=2-> (    -8.960,      9.500)
!                       opt=3-> ( 9.140,  0.000)
!                      opt=4-> ( 0.5946107375736628, 0.6855384522353541)
!                  flag='integer: 1-> no info flag, 0 or nothing-> info flag'
!          )
    subroutine mat_dump_vec(vec,text,opt,flag)
        implicit none
        !
        class(*), dimension(:), intent(in) :: vec
        character(*), intent(in), optional :: text
        integer, intent(in), optional :: opt
        integer, intent(in), optional :: flag
        !system
        integer*4 N
        character(50) fmtStringCmplx,fmtString,fmtStringInt,typeof
        !
        if (present(opt))then!può essere aumentato in futuro
           select case (opt)
                case(1)
                    fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
                    fmtString = '(9999("(",ES12.4,")",:,1X))'
                case(2)
                    fmtStringCmplx = '(9999("(",F10.3,",",1X,F10.3,")",:,1X))'
                    fmtString = '(9999("(",F10.3,")",:,1X))'
                case(3)
                    fmtStringCmplx = '(9999("(",F3.0,",",1X,F3.0,")",:,1X))'                      
                    fmtString = '(9999("(",F20.10,")",:,1X))'
                case(4)
                    fmtStringCmplx = '(9999("(",F20.16,",",1X,F20.16,")",:,1X))'                      
                    fmtString = '(9999("(",F20.16,")",:,1X))'
                case(5)!mathematica
                    fmtStringCmplx = '("{",9999(F20.16,"+ I(",1X,F20.16,"),",:,1X))'                      
                    fmtString = '("{",9999("(",F20.16,")",:,1X),"}")'
                case default
                    fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
                    fmtString = '(9999("(",ES12.4,")",:,1X))'
           endselect
        else
            fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
            fmtString = '(9999("(",ES12.4,")",:,1X))'
        endif
        fmtStringInt = '(9999("(",I0,")",:,1X))'    
        !print*,''
        print*,text
        typeof='unknown'
        !
        N=size(vec)
        select type(vec)
            type is (complex(16))
                typeof='complex*16'
                write(*,fmt=fmtStringCmplx) vec
            type is (complex(8))
                typeof='complex*8'
                write(*,fmt=fmtStringCmplx) vec
            type is (complex(4))
                typeof='complex*4'
                write(*,fmt=fmtStringCmplx) vec
            type is (integer(2))
                typeof='integer*2'
                write(*,fmt=fmtStringInt) vec
            type is (integer(4))
                typeof='integer*4'
                write(*,fmt=fmtStringInt) vec
            type is (real(4))
                typeof='real*4'
                write(*,fmt=fmtString) vec
            type is (real(8))
                typeof='real*8'
                write(*,fmt=fmtString) vec
            type is (real(16))
                typeof='real*16'
                write(*,fmt=fmtString) vec
            type is (character(len=*))
                typeof='character'
                write(*,'(A)') vec
            type is (logical)
                typeof='logical'
                write(*,fmt='(9999("(",L,")",:,1X))') vec       
        endselect                       
        if(present(flag))then
            select case (flag)
            !no info message
                case(1)                   
                    print*,''
                case default
                    print*,'attenzione: flag/=1'
            endselect
        else
            write(*,'(" [vector details: ",A,", dim: ",I0,"]")') trim(typeof),N
        endif
        !print*,''
        return
    endsubroutine mat_dump_vec

    subroutine mat_dump_vec_2_file(vec,unit_file,text,opt,flag)
        implicit none
        !
        class(*), dimension(:), intent(in) :: vec
        integer, intent(in) :: unit_file
        character(*), intent(in), optional :: text
        integer, intent(in), optional :: opt
        integer, intent(in), optional :: flag
        !system
        integer*4 N
        character(50) fmtStringCmplx,fmtString,fmtStringInt,typeof
        !
        if (present(opt))then!può essere aumentato in futuro
           select case (opt)
                case(1)
                    fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
                    fmtString = '(9999("(",ES12.4,")",:,1X))'
                case(2)
                    fmtStringCmplx = '(9999("(",F10.3,",",1X,F10.3,")",:,1X))'
                    fmtString = '(9999("(",F10.3,")",:,1X))'
                case(3)
                    fmtStringCmplx = '(9999("(",F3.0,",",1X,F3.0,")",:,1X))'                      
                    fmtString = '(9999("(",F20.10,")",:,1X))'
                case(4)
                    fmtStringCmplx = '(9999("(",F20.16,",",1X,F20.16,")",:,1X))'                      
                    fmtString = '(9999("(",F20.16,")",:,1X))'
                case(5)!mathematica
                    fmtStringCmplx = '("{",9999(F20.16,"+ I(",1X,F20.16,"),",:,1X))'                      
                    fmtString = '("{",9999("(",F20.16,")",:,1X),"}")'
                case default
                    fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
                    fmtString = '(9999("(",ES12.4,")",:,1X))'
           endselect
        else
            fmtStringCmplx = '(9999("(",ES12.4,",",1X,ES12.4,")",:,1X))'
            fmtString = '(9999("(",ES12.4,")",:,1X))'
        endif
        fmtStringInt = '(9999("(",I0,")",:,1X))'    
        !print*,''
        write(unit_file,*) text
        typeof='unknown'
        !
        N=size(vec)
        select type(vec)
            type is (complex(16))
                typeof='complex*16'
                write(unit_file,fmt=fmtStringCmplx) vec
            type is (complex(8))
                typeof='complex*8'
                write(unit_file,fmt=fmtStringCmplx) vec
            type is (complex(4))
                typeof='complex*4'
                write(unit_file,fmt=fmtStringCmplx) vec
            type is (integer(2))
                typeof='integer*2'
                write(unit_file,fmt=fmtStringInt) vec
            type is (integer(4))
                typeof='integer*4'
                write(unit_file,fmt=fmtStringInt) vec
            type is (real(4))
                typeof='real*4'
                write(unit_file,fmt=fmtString) vec
            type is (real(8))
                typeof='real*8'
                write(unit_file,fmt=fmtString) vec
            type is (real(16))
                typeof='real*16'
                write(unit_file,fmt=fmtString) vec
            type is (character(len=*))
                typeof='character'
                write(unit_file,'(A)') vec
            type is (logical)
                typeof='logical'
                write(unit_file,fmt='(9999("(",L,")",:,1X))') vec       
        endselect                       
        if(present(flag))then
            select case (flag)
            !no info message
                case(1)                   
                    print*,''
                case default
                    print*,'attenzione: flag/=1'
            endselect
        else
            write(unit_file,'(" [vector details: ",A,", dim: ",I0,"]")') trim(typeof),N
        endif
        !print*,''
        return
    endsubroutine mat_dump_vec_2_file
!LengthString(
!              STR(character(*)) = 'string one wants to fine the length'
!          )
!fined the length of a string, .e.g. for a char(3) STR = 'a' returns 1
    FUNCTION LengthString(STR)
        IMPLICIT NONE
        CHARACTER(*), INTENT(IN) :: STR
        INTEGER(8)               I, LengthString
        
        LengthString = -1 !default is -1 just because.
        
        DO I=LEN(STR), 1, -1
        IF(STR(I:I) .ne. ' ') THEN
            LengthString = I
            EXIT
        ENDIF
        ENDDO

        RETURN
    ENDFUNCTION LengthString

!split(
!              input_line,
!              array  =  output array
!              delimiters = delimiters
!          )
!returns the array with the elements separeted by delimeter
    subroutine split(input_line,array,delimiters,order,nulls)
        !! given a line of structure " par1 par2 par3 ... parn " store each par(n) into a separate variable in array.
        !!
        !! * by default adjacent delimiters in the input string do not create an empty string in the output array
        !! * no quoting of delimiters is supported
        character(len=*),intent(in)              :: input_line  !! input string to tokenize
        character(len=*),optional,intent(in)     :: delimiters  !! list of delimiter characters
        character(len=*),optional,intent(in)     :: order       !! order of output array sequential|[reverse|right]
        character(len=*),optional,intent(in)     :: nulls       !! return strings composed of delimiters or not ignore|return|ignoreend
        character(len=:),allocatable,intent(out) :: array(:)    !! output array of tokens
    
        integer                       :: n                      ! max number of strings INPUT_LINE could split into if all delimiter
        integer,allocatable           :: ibegin(:)              ! positions in input string where tokens start
        integer,allocatable           :: iterm(:)               ! positions in input string where tokens end
        character(len=:),allocatable  :: dlim                   ! string containing delimiter characters
        character(len=:),allocatable  :: ordr                   ! string containing order keyword
        character(len=:),allocatable  :: nlls                   ! string containing nulls keyword
        integer                       :: ii,iiii                ! loop parameters used to control print order
        integer                       :: icount                 ! number of tokens found
        integer                       :: ilen                   ! length of input string with trailing spaces trimmed
        integer                       :: i10,i20,i30            ! loop counters
        integer                       :: icol                   ! pointer into input string as it is being parsed
        integer                       :: idlim                  ! number of delimiter characters
        integer                       :: ifound                 ! where next delimiter character is found in remaining input string data
        integer                       :: inotnull               ! count strings not composed of delimiters
        integer                       :: ireturn                ! number of tokens returned
        integer                       :: imax                   ! length of longest token
    
        ! decide on value for optional DELIMITERS parameter
        if (present(delimiters)) then                                     ! optional delimiter list was present
            if(delimiters.ne.'')then                                       ! if DELIMITERS was specified and not null use it
                dlim=delimiters
            else                                                           ! DELIMITERS was specified on call as empty string
                dlim=' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0) ! use default delimiter when not specified
            endif
        else                                                              ! no delimiter value was specified
            dlim=' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0)    ! use default delimiter when not specified
        endif
        idlim=len(dlim)                                                   ! dlim a lot of blanks on some machines if dlim is a big string
    
        if(present(order))then; ordr=order; else; ordr='sequential'; endif ! decide on value for optional ORDER parameter
        if(present(nulls))then; nlls=nulls; else; nlls='ignore'    ; endif ! optional parameter
    
        n=len(input_line)+1                        ! max number of strings INPUT_LINE could split into if all delimiter
        allocate(ibegin(n))                        ! allocate enough space to hold starting location of tokens if string all tokens
        allocate(iterm(n))                         ! allocate enough space to hold ending location of tokens if string all tokens
        ibegin(:)=1
        iterm(:)=1
    
        ilen=len(input_line)                                           ! ILEN is the column position of the last non-blank character
        icount=0                                                       ! how many tokens found
        inotnull=0                                                     ! how many tokens found not composed of delimiters
        imax=0                                                         ! length of longest token found
    
        select case (ilen)
    
        case (0)                                                      ! command was totally blank
    
        case default                                                   ! there is at least one non-delimiter in INPUT_LINE if get here
            icol=1                                                      ! initialize pointer into input line
            INFINITE: do i30=1,ilen,1                                   ! store into each array element
                ibegin(i30)=icol                                         ! assume start new token on the character
                if(index(dlim(1:idlim),input_line(icol:icol)).eq.0)then  ! if current character is not a delimiter
                iterm(i30)=ilen                                       ! initially assume no more tokens
                do i10=1,idlim                                        ! search for next delimiter
                    ifound=index(input_line(ibegin(i30):ilen),dlim(i10:i10))
                    IF(ifound.gt.0)then
                        iterm(i30)=min(iterm(i30),ifound+ibegin(i30)-2)
                    endif
                enddo
                icol=iterm(i30)+2                                     ! next place to look as found end of this token
                inotnull=inotnull+1                                   ! increment count of number of tokens not composed of delimiters
                else                                                     ! character is a delimiter for a null string
                iterm(i30)=icol-1                                     ! record assumed end of string. Will be less than beginning
                icol=icol+1                                           ! advance pointer into input string
                endif
                imax=max(imax,iterm(i30)-ibegin(i30)+1)
                icount=i30                                               ! increment count of number of tokens found
                if(icol.gt.ilen)then                                     ! no text left
                exit INFINITE
                endif
            enddo INFINITE
    
        end select
    
        select case (trim(adjustl(nlls)))
        case ('ignore','','ignoreend')
            ireturn=inotnull
        case default
            ireturn=icount
        end select
        allocate(character(len=imax) :: array(ireturn))                ! allocate the array to return
        !allocate(array(ireturn))                                       ! allocate the array to turn
    
        select case (trim(adjustl(ordr)))                              ! decide which order to store tokens
        case ('reverse','right') ; ii=ireturn ; iiii=-1                ! last to first
        case default             ; ii=1       ; iiii=1                 ! first to last
        end select
    
        do i20=1,icount                                                ! fill the array with the tokens that were found
            if(iterm(i20).lt.ibegin(i20))then
                select case (trim(adjustl(nlls)))
                case ('ignore','','ignoreend')
                case default
                array(ii)=' '
                ii=ii+iiii
                end select
            else
                array(ii)=input_line(ibegin(i20):iterm(i20))
                ii=ii+iiii
            endif
        enddo
    endsubroutine split
!strip_spaces(
!              string = input string,
!              outstring = output string,
!          )
!return the string without blank spaces
    subroutine strip_spaces(string, outstring)
        character(len=*)                             :: string
        character(len=:), allocatable, intent(inout) :: outstring
        integer                                      :: stringLen 
        integer                                      :: last, actual
    
        stringLen = len (string)
        last = 1
        actual = 1
    
        do while (actual < stringLen)
            if (string(last:last) == ' ') then
                actual = actual + 1
                string(last:last) = string(actual:actual)
                string(actual:actual) = ''
            else
                last = last + 1
                if (actual < last) &
                    actual = last
            endif
        end do
        allocate(character(len_trim(string)) :: outstring)
        outstring = trim(string)
        
        ! allocate(outstring(len_trim(string)))
    
    endsubroutine strip_spaces
endmodule debugging