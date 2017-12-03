function [Ch_Rf_Des,Ch_Rf,Ch_Rf_a]=Gen_Channel(mode,Nbr_usager,Long_Message,Long_Slot,PG,v)

		% Crée le 07/10/2010 par B.Belgacem

%/------------------- Début fonction -----------------------\%


switch mode

    case {1}

         Ch_Rf_Des   = ones(Nbr_usager,Long_Message*PG) ;
 
        Ch_Rf       = ones(Nbr_usager,Long_Message) ;
 
        Ch_Rf_a     = ones(Nbr_usager,Long_Message/Long_Slot );

    case {2}

         Ch_Rf_Des   = zeros(Nbr_usager,Long_Message*PG) ;
 
        Ch_Rf       = zeros(Nbr_usager,Long_Message) ;
 
        Ch_Rf_a     = zeros(Nbr_usager,Long_Message/Long_Slot );

        %--------------------------------
        
for i=1:Nbr_usager       
	Ch_Rf(i,:) = Gen_Rayleigh_Ch(Long_Message,PG,v(i));

end
y=1;
        
for j=1:Long_Message,

            for p=1:PG,
                Ch_Rf_Des(:,y) = Ch_Rf(:,j);

                y=y+1;

            end


end

y=1;
for i=1:Long_Message/Long_Slot

            for j=1:Long_Slot

                Ch_Rf_a(:,i) = Ch_Rf_a(:,i)+Ch_Rf(:,y);

                y=y+1;

            end
        
end

        Ch_Rf_a = 1/Long_Slot*Ch_Rf_a;

end;

