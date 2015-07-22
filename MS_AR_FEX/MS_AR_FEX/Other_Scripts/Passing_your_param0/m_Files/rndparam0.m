function param0=rndparam0(x,ar,k)
    
    for ik=1:k
        Spec.AR0(:,ik)=randn(1,ar);  % param0 for Ar parameters with normal rnd (mean=0, std=1)
        Spec.const0(1,ik)=mean(x)+randn(1,1)*std(x);  % param0 for mean constant at state 1 with mean=mean(x) and std=std(x)
        Spec.Std0(1,ik)=rand(1,1);  % param0 for model's standard deviation at state 1
    end
    
    R=rand(k);
    
    Spec.p0=R./repmat(sum(R),k,1); % Making sure the Spec.p0 sums 1 for each collum
    
    param0=spec2param0(Spec,ar,k); % converting Spec to param0 vector
    
    


