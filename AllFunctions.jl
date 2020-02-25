module AllFunctions

function div_n_ab( n::Int, lo::Int64 = 1, hi::Int64 = n )                       
#=
    Divisores del numero n entre los rangos "lo" and "hi"
=#
    ρ = collect( 1:Int( floor( sqrt( n ) ) ) ); # los numeros de 1 en 1 de la raiz cuadrada 
    σ1 = findall( n.%ρ .== 0 ); # divisores de la raiz cuadrada ( residuo = 0 )
    σ2 = Int.( ( n )./( σ1 ) ); # Sacar los pares ( de 100, 2-50, 10-10, etc.)
    σ = sort( unique( vcat( σ1, σ2 ) ) ); # remover duplicados, concatenar, ordenar
    aux1 = @isdefined lo;
    aux2 = @isdefined hi;
    if aux1 && aux2
        rn = σ[ findall( hi .>= σ .>= lo ) ];
        return rn
    else
        return σ
    end
end

function divs0prime( number::Int )
    divs = AllFunctions.div_n_ab(number)
    divs = filter!( e -> e ≠ 1, divs ); # removiendo los limites
    divs = filter!( e -> e ≠ number, divs ); # removiendo los limites
    return divs
end

function OS( OS::String )
# para los paths de mi lap o del lab
    if OS == "linux"
        com = "/"
    elseif OS == "windows"
        com = "\\"
    end
    return com
end

function checkpath( workpath::String )
# checar si el folder existe, si no, hacerlo
    if isdir( workpath ) == false
        mkpath( workpath )
    end
end

function find_files( path_main::String, key::String, Γ::String ) # depende de OS!!!
# buscar en el directorio Path_main, los archivos que terminen con key
    searchdir(path_main::String, key::String) = filter( x -> endswith(x, key), readdir(path_main) );
    files = searchdir( path_main, key );
    if path_main[ end ] == collect( Γ )[ 1 ]
        files = string.( path_main, files );
    else
        files = string.( path_main, Γ, files );
    end
    return files
end

function find_dirs( path_main::String, Γ::String ) # depende de OS!!!
# buscar en el directorio Path_main, los subdirectorios existentes
    dirs = filter(x -> isdir(joinpath(path_main, x)), readdir(path_main))
    if path_main[ end ] == collect( Γ )[ 1 ]
        dirs = string.( path_main, dirs );
    else
        dirs = string.( path_main, Γ, dirs );
    end
    return dirs
end

function neighborgs( center::Int, d::Int )
# obtienen la "d"-vecindad del canal "center" #
    A = reshape( 1:4096, 64, 64 );
    A = Int.( A' );
    x_c = findall( A .== center )[ ][ 2 ]
    y_c = findall( A .== center )[ ][ 1 ]
    aux = [ ( x_c - d ),( x_c + d ), ( y_c - d ), ( y_c + d ) ]
    aux[ aux .< 1 ] .= 1;
    aux[ aux .> 64 ] .= 64;
    neigh = A[ aux[ 3 ]:aux[ 4 ], aux[ 1 ]:aux[ 2 ] ];
    return neigh
end

function A_minus_B( A::Array, a::Array )
    #=
    Compara cada fila del array chico a con las filas del array grande A
    y si están repetidas las elimina del array grande A
    Devuelve el array grande sin las filas del chico
    =#
    n = size( A, 2 );
    m = size( a, 2 );
    if n == m
        for i = 1:size( a , 1 )
            B = reshape( a[ i, : ],( 1, m ) );
            C = ( A .== B );
            C = sum( C, dims = 2 ) .== size( C, 2 );
            D = collect( 1:size( C, 1 ) );
            co = filter!( e -> e .!= 0, vec( Int.( D.*( 1 .- C ) ) ) ) ;
            A = A[ co, : ];
        end
    else
        println("array chico no tiene mismo numero de columnas que array grande")
    end
    return A
end

function sats(data::Array, HIthr::Int, LOthr::Int)
    saturados = vcat(findall(data .>= HIthr),findall(data .<= LOthr));
    ChFrSat = zeros( Int, size( saturados, 1 ), 2 ); # preallocation
    if !isempty( saturados )
        for j = 1:size( saturados, 1 )
            ChFrSat[ j, 1 ] = saturados[ j ].I[ 1 ]; # channel
            ChFrSat[ j, 2 ] = saturados[ j ].I[ 2 ]; # Frame
        end
    end
    PSP = round( 
        ( size( ChFrSat, 1 )*100 )/( size( data, 1 )*size( data, 2 ) ), digits = 2 
        );
    println( string( " Hay ", PSP,"% de saturación." ) );
    return ChFrSat
end

function ESU( ch::Array, thr::Real, d::Int )
    #= Se obtienen los frames y voltajes que sobrepasan el umbral establecido. Si hay eventos 
    supraumbral a d-frames de distancia, se selecciona aquel que tenga menor voltaje.
    =#
    OK = findall( ch .<= thr );
    if !isempty(OK)
        init = 1;
        while init == 1
            A = OK[ 1:( end - 1 ) ]; B = OK[ 2:end ];
            C = B .- A;
            D = findall( C .<= d ) .+ 1;
            if isempty( D )
                init = 0;
            else
                nook = zeros( Int, size( D, 1 ) );
                for E = 1:size( D, 1 )
                    # si el primero es menor que el segundo, quita el segundo
                    if isless( (ch[ OK[ D[ E ] ] ]), (ch[ OK[ D[ E ] - 1 ] ] ) ) 
                        nook[ E ] = D[ E ] - 1;
                    else
                        nook[ E ] = D[ E ];
                    end
                end
                if size( nook, 1 ) > 1
                    OK[ unique( nook ) ] .= 0;
                    filter!( x -> x != 0, OK );
                else
                    OK = OK[ Bool.( OK .!= OK[ nook[ 1 ] ] ) ];
                end
            end
        end
        channelC = zeros( size( OK, 1 ), 2 );
        channelC[ :, 1 ] = OK;
        channelC[ :, 2 ] = ch[ OK ];
    else
        channelC = [ ];
    end
    return channelC
end

# -------------------------------------------------------------------------------------------- !!!
export div_n_ab
export OS
export checkpath 
export find_files 
export neighborgs 
export find_dirs 
export divs0prime
export A_minus_B
export sats
export ESU

end
                           
                           
                           
                           
