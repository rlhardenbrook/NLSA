function tag = getDefaultTag( obj )
% GETDEFAULTTAG  Get default tag of nlsaKoopmanOperator_diff objects
%
% Modified 2020/08/01

tag = getDefaultTag@nlsaKoopmanOperator( obj );

tag = [ tag sprintf( 'diff_%s_eps%1.3g', getRegularizationType( obj ), ...
                                         getRegularizationParameter( obj ) ) ];
