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

using HDF5
function brw_things( file_brw::String ) # depende de brw file
	brw = h5open( file_brw, "r" );
	# channels and coords for each one
	Chs = read( brw[ "/3BRecInfo/3BMeaStreams/RawRanges" ] )[ "Chs" ]; 
	x = zeros( Int, size( Chs, 1 ), 3 );
	for i = 1:size( Chs, 1 )
	   # Channels, [number coord1 coord2]
	   x[ i, : ] = [ i, Chs[ i ].data[ 1 ], Chs[ i ].data[ 2 ] ]; 
	end  
	VC = read( brw[ "/3BData/3BInfo/3BNoise/ValidChs" ] ); # Valid Channels (software desition)
	y = zeros( Int, size( VC, 1 ), 3 );
	for i = 1:size( VC, 1 )
	# Valid Channels, [number coord1 coord2]
	   y[i,:] = [ 
	   ( ( VC[ i ].data[ 1 ] - 1)*64 + VC[ i ].data[ 2 ] ), VC[ i ].data[ 1 ], VC[ i ].data[ 2 ] 
	   ];
	end
	# numero real de frames registrados
	RawEncodedTOC = read( brw[ "/3BData/RawEncodedTOC" ] )[ 2 ]; 
	RecVars = read( brw[ "/3BRecInfo/3BRecVars" ] );
	SI = RecVars[ "SignalInversion" ][ ];
	mV = RecVars[ "MinVolt" ][ ];
	MV = RecVars[ "MaxVolt" ][ ];
	BD = RecVars[ "BitDepth" ][ ];
    Offset = SI*mV;
    Factor = SI*( MV - mV )/( 2^BD );
    NRecFrames = RecVars[ "NRecFrames" ][ ]; # numero propuesto de cuadros registrados
    SamplingRate = RecVars[ "SamplingRate" ][ ];
    dset = brw[ "3BData/RawEncoded" ];
  
	vars = Dict(
        "Offset"        => Offset,
        "Factor"        => Factor,
        "NRecFrames"    => NRecFrames,
        "SamplingRate"  => SamplingRate,
        "Chs"   	    => x,
		"VC" 	        => y,
        "RawEncodedTOC" => RawEncodedTOC,
        "dset"          => dset,
        "MaxVolt"       => MV,
        "MinVolt"       => mV
        
    );
    return vars
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
# -------------------------------------------------------------------------------------------- !!!
export div_n_ab, OS, checkpath, brw_things, find_files, neighborgs, find_dirs, divs0prime, A_minus_B, sats

end
                           
                           
                           
                           
