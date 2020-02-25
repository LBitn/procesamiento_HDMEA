module BRW_things

using HDF5
function brw_things( file_brw::String ) # depende de brw file
	brw = h5open( file_brw, "r" );
	# channels and coords for each one
	Chs = read( brw[ "/3BRecInfo/3BMeaStreams/RawRanges" ] )[ "Chs" ]; 
    
    x = zeros( Int, size( Chs, 1 ), 3 );
	for i = 1:size( Chs, 1 )
	   x[ i, : ] = [ i, Chs[ i ].data[ 1 ], Chs[ i ].data[ 2 ] ]; # Channels, [number coord1 coord2]
	end
 
    dset = brw[ "3BData/RawEncoded" ];
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
    BRWname = dset.file.filename;
    
	vars = Dict(
        "Offset"        => Offset,
        "Factor"        => Factor,
        "NRecFrames"    => NRecFrames,
        "SamplingRate"  => SamplingRate,
        "Chs"   	    => x,
        "RawEncodedTOC" => RawEncodedTOC,
        "BRWname"       => file_brw,
        "MaxVolt"       => MV,
        "MinVolt"       => mV,
        "dset"          => dset
    );
    return vars
end
export brw_things
end
