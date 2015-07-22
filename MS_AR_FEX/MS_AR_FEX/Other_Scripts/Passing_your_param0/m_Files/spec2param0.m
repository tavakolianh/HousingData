function [param0]=spec2param0(Spec,ar,k)
    
    if nargin()~=3
        error('there should be 3 arguments at spec2param0()');
    end
    
    if any( (sum(Spec.p0)>1.0001)|(sum(Spec.p0)<.999) )
        error('The sum of each collum in Spec.p0 should be equal to 1 (they are probabilities..)')
    end
    
    if any(Spec.p0<0+Spec.p0>1)
        error('The Spec.p stores probabilities and they should be lower than 1 and higher than 0, unless you have a new theory about probabilities (it should be a very convincing one)');
    end
    
    if size(Spec.AR0,2)~=k
        error('The Spec.AR0 has a different number of collum than the value of k.');
    end
    
    if size(Spec.Std0,2)~=k
        error('The Spec.Std has a different number of collum than the value of k.');
    end

    if size(Spec.const0,2)~=k
        error('The Spec.const0 has a different number of collum than the value of k.');
    end
    
    if any(Spec.Std0<0)
        error('The Spec.Std0 is a standard deviation and should be positive. Ohhhhh wait, are you saying that you''ve created the most perfect model, where the error is replaced by a varying degree of certanty ? Thats nice, you should try publish it...');
    end
    
    param0=zeros(1,k+k+ar*k+k^2);
    
    param0(1:k)=Spec.Std0(1,1:k);
    
    param0(k+1:2*k)=Spec.const0(1,1:k);

    for i=0:k-1
        param0(2*k+1+i*ar:k+k+ar+i*ar)=Spec.AR0(1:ar,i+1);
        param0(end-k^2+1+i*k:end-k^2+k+i*k)=Spec.p0(i+1,1:k);
    end  