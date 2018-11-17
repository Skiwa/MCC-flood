/**
* Name: basic
* Author: carole
* Basic flood model for WIC project, adapted from SPRITE model
*
* Description: 
* Tags: Tag1, Tag2, TagN
*/

model basic


global
{
	/* ************************************************************************ *
	 * **************************** VARIABLES AGENT MONDE ********************* *
	 * ************************************************************************ */
	  
	//dossier contenant les fichiers à lire
	string dossierFichiers <-  "../includes/" ;
	string dossierImages <-  "../images/" ;
		
	//fichier à lire (grille Oleron) ou en écriture en god mode.
	string cell_grid <- dossierFichiers +"grille_oleron.csv";
	string cell_grid_mod <- dossierFichiers +"grille_oleronP.csv";
	
	//image des boutons
	list<image_file> images <- [		
		image_file(dossierImages +"titre0.png"),
		image_file(dossierImages +"titre1.png"),
		image_file(dossierImages +"titre2.png"),
		image_file(dossierImages +"impots.png"),
		image_file(dossierImages +"PLU_Urb.png"),
		image_file(dossierImages +"BA.png"),
		image_file(dossierImages +"promotion.png"),
		image_file(dossierImages +"PLU_Agri.png"),
		image_file(dossierImages +"rocher.png"),
		image_file(dossierImages +"produits.png"),
		image_file(dossierImages +"PLU_nat.png"),
		image_file(dossierImages +"dune.png"),
		image_file(dossierImages +"dommage.png"),
		image_file(dossierImages +"exproprier.png"),
		image_file(dossierImages +"reparer.png"),
		image_file(dossierImages +"affichage.png"),
		image_file(dossierImages +"batiment.png"),
		image_file(dossierImages +"reparer_all.png"),
		image_file(dossierImages +"no_action.png"),
		image_file(dossierImages +"blanc.png"),
		image_file(dossierImages +"demolir.png")
	]; 


	//chargement du contour de l'ile et du cadre de mer
	geometry shape <- envelope(file(dossierFichiers+"rect.shp"));

	
	//******quelques définitions de listes permettant de manipuler des parcelles précises*******
	list<parcelle> sea_cells;				//parcelles de mer
   	list<parcelle> cells_close_to_sea;		//parcelles proches de la mer
	list<parcelle> earth_cells;				//liste des parcelles de terre
	//****************************************


	//***************** parametres de submersion*****************
	//définit pour chaque submersion sa force (de 1 pas fort à 10 très fort)
	int submersion_intensity <- rnd(10);
			
	//string weather_info<-"";
	//coefficient de maree
	int tidal_coeff<-0;
	//vitesse vent						
	float wind_speed;
	string wind_info <- "";
	string wave_info <- "";
	string weather_info <- "";
	string water_info <- "";
	
	//pour choisir l'orientation du vent dans les scénarios
	int submersion_ww <- rnd(0,7);	
	int direction <- rnd(1,8);			// direction entre 1 et 8, parameter
	int wind_row<-0;						//colonne pour direction du vent
	int wind_line<-0;					//ligne pour direction du vent
	float atmosph_pressure<-0.0;			//pression atmosphérique
	float added_water<-0.0; 				//eau en plus sur la mer lors d'une submersion
	float normal_water_heigth <- 5.2;  	//hauteur normal de l'eau en marrée haute (sert de référence)
	float coeff_wave <-0.3;   			//rapport entre vent en beaufort et hauteur de vague (ex:si 10 Bft et coeff_wave=0.5 => vague de 5m)
	
	bool end_flowing<-false;				//annonce fin de la submersion (repasse en mode tour normal)
   	list<parcelle> cells_to_flow;			//parcelles avec de l'eau et de la terre a coté
	list<string> storm_name <- ["Carole", "Odile", "Mira", "Germaine", "Raymonde", "Corinne", "Gertrude", "Cunégonde", "Martine", "Jacqueline"];
	int sub_nb<-rnd(0,9); 					//iteration n° de submersion
	int sub_year update: 2010+sub_nb;
	
	//int last_sub_year<-2010;										--> sub_year
	//string last_sub_name<-"Xynthia";									--> storm_name[sub_nb]
	//string last_sub_waether_P<-"Pression atmosphèrique : 969 hPa"; 			--> atmosph_pressure
	//string last_sub_waether_M<-"Coefficient de marée : 102"; 				--> tidal_coeff
	//string last_sub_waether_V<-"Vent de Sud Ouest d'une force de 10 Bf";		--> weather_info
	
	//float total_dom;
	//float total_dam_hist<-0.0;
	

	

	/*******************************************************************		
	 ****************** INITIALISATION AGENT MONDE *********************
	 *******************************************************************/
	init 	{			
		//création des presidents 
		do init_cells;														//initialisation des parcelles de la grille a partir de la matrice
		//SIMPLIF// do init_butons;
	}
	//**************************************
	// 				FIN DU INIT
	//************************************
	

		
	//initialisation des cellules de la grille a partir du csv d'Oléron
	action init_cells
	{
		//passage par un agent cell qui réceptionne les données et les renvoie 
		matrix init_data <- matrix(csv_file(cell_grid));
		
		int pop <-1;		
		//pop est un incrément récupérant les n° de colonne dans les donnees
		//on copie les infos de chaque parcelle à partir de la matrice
		ask parcelle
		{
			pop <- init_parcel(init_data,pop);
		}
		
		// AFTER all initialised
		ask parcelle {
			//sera utilisé pour la submersion
			earth_neighbour_cells <- neighbour_cells where (!each.is_water);
			is_earth_neighbour_cells<- (length(earth_neighbour_cells)>0);
		}
		
		// global variables	related to sea
		sea_cells <- parcelle where each.is_sea;
		earth_cells <- parcelle - sea_cells;
		
		// AFTER sea parcels initialised
		ask parcelle {	
			if !is_sea {
				sea_neighbour_cells<- neighbour_cells where (each.is_sea);  //remplissage liste de parcelles de mer à coté et 
				close_sea <- (length(sea_neighbour_cells)>0);			// remember the sea is close if at least one neighbour cell is sea	
				closestSea <- (sea_cells closest_to(self)); 
				using topology(world) {	distance_to_sea <- self distance_to closestSea;}		//Calcul de la distance en m de la parcelle à la mer (reste inchangé pendant toute la partie)
			}
		}
	} //end action	


	
	
	//*************************************************
	//			DYNAMIQUE DU MODELE
	//*************************************************
	
	// simplifie a une seule des DEUX dynamiques:
	// * si NT=false:phase de submersion (avec sa propre dynamique qui s'arrete quand end_flowing)
	
	// Action pour déclencher la submersion - juste pour debugge ; à ne pas inclure dans le jeu

	// doit etre un reflexe pour se faire automatiquement au debut de chaque tour
	reflex submersion_reflex {
		if (end_flowing){ 
			do restore_world;
			end_flowing<-false;
			do pause;
		}
		else {
			do submerge;
		}	
	}

	reflex submersion when: cycle=1 {
		do init_sub;
	}

	user_command "Force submersion !" {
		do init_sub;	
	}

	user_command "Restart" {
		do restore_world;
	}
	
	user_command "Sauvegarder carte" {
		save ["altitude","is_sea","digue","densite_bati","obstacle_height","constructible","historique"] rewrite:false to:cell_grid_mod type:"csv";			//save des indicateurs en début de chaque année	
		ask parcelle {
			save [altitude,is_sea,dyke,dyke_height,water_height_history] rewrite:false to:cell_grid_mod type:"csv";	
		}
		write ("Carte sauvegardée sous grille_oleron2.csv");
	}	
	
	//action qui gère la submersion
	action submerge {		
		//on choisit les parcelles à faire déverser (celles avec de l'eau, à coté d'une cellule de terre, et ne s'étant pas encore déversée)
		cells_to_flow <- parcelle where(each.is_water and each.is_earth_neighbour_cells and !each.already);

		//write ("Submersion, cells to flow = "+cells_to_flow);

		//on demande à ces cellules de se derverser
	  	ask cells_to_flow 			{	
			do flow;				//déversion (hors vague)
			row <- grid_x;		//utilisé pour les vagues
			line <-grid_y;		//utilisé pour les vagues
			do wave_wind;		//propagagtion de vagues
			earth_neighbour_cells <- neighbour_cells where (!each.is_water);
		}
		
		// condition de fin (s'il n'y a plus de parcelle pour le derversement, ça s'arrete)
		if length(cells_to_flow)=0 	{	
			write("end of submersion at cycle = "+cycle);
			end_flowing<-true;
		}
	} //fin de l'action submerge
	
	
	// RAZ des parcelles après une submersion
	action restore_world {
		write("Restoring world at cycle "+cycle);
		//pour les parcelles de terre
		ask (earth_cells) {	 	
			is_water<-false;
			already<-false;
			water_height<-0.0;
			water_altitude <-altitude;
		}
		
		//pour les parcelles de mer
		ask (sea_cells){
			water_height<-(-altitude);
			water_altitude <-0.0;
			already<-false;
			is_water<-true;
		}
		
		//réinitialisation de ces listes pour la prochaine submersion
		ask parcelle {
			earth_neighbour_cells <- neighbour_cells where (!each.is_water);
			is_earth_neighbour_cells <- (length(earth_neighbour_cells)>0);
		}
	}

	

	// action qui oriente le sens de propagation des vagues en fonction de la direction du vent
	action init_wind
	{
		write("global init wind");
		
		//string dir<-"";
		//direction<-7;  //Xynthia
		list<int> wind_rows <-  [1,  0, -1, 0,  1, -1, -1, 1];
		list<int> wind_lines <- [0, -1,  0, 1, -1, -1,  1, 1];
		list directions <- ["d'Est", "de Nord", "d'Ouest", "de Sud", "de Nord Est", "de Nord Ouest", "de Sud Ouest", "de Sud Est"] ;
		//impact de la direction (index entre 1 et 8)
		wind_row <- wind_rows[direction];
		wind_line<-wind_lines[direction];
		//dir<-directions[direction];			
		wind_info<-"Vent "+directions[direction];
		
		//int submersion_intens <- rnd(10);
		int WS <- int(min([12,max([0,6+submersion_intensity/2])]));  //force du vent en beaufort
		//WS<-10;		//Xynthia
		float WH <- WS*coeff_wave; //hauteur de vague
		wind_info<-wind_info+" d'une force de "+WS+" Beaufort\n";
		wave_info <- "Hauteur de vague:"+WH with_precision 2+" m\n";
		//last_sub_waether_V<-"Vent de "+dir+" d'une force de "+WS+" Bf";
		//initialisation des vagues en mer
		ask sea_cells	
		{
			wave_height<-WH;
			wave_altitude <-water_altitude+wave_height;
		}
		write("wind info after init = "+wind_info);
	}
	

	//action qui définit la hauteur d'eau de submersion en fonction de la pression atmosphérique et du coefficient de marée
	action init_sub
	{
		write("global init sub");
		string storm_name_year <-storm_name[sub_nb];
		weather_info<-"Informations météo sur la tempête "+storm_name_year+":\n";
		atmosph_pressure<-1015-(submersion_intensity*4.0);  //pression tirée au hasard dans une gamme plausible (pour une tempete)
		//atmosph_pressure<-969;   //Xynthia
		tidal_coeff<-90+submersion_intensity*2;    			//coeff tiré au hasard correspondant à de grande marées
		//tidal_coeff<-102;  //Xynthia
		added_water<-0.01*(1015-atmosph_pressure)+(1.2+tidal_coeff/100)*2.5-normal_water_heigth;		//formule issue de la littérature
		ask sea_cells	{ water_altitude <- added_water;	}
		weather_info <- weather_info+
						"Pression atmosphèrique:"+atmosph_pressure+" hPa \n"+
						"Coefficient de marée:"+tidal_coeff+"\n"+
						"Hauteur d'eau:"+added_water with_precision 2+" m\n";
		water_info <- "Hauteur d'eau:"+added_water with_precision 2+" m\n";
		do init_wind;
		sub_nb<-(sub_nb+1);
		
		//last_sub_name<-storm_name_year;
		//last_sub_waether_P<-"Pression atmosphèrique:"+atmosph_pressure+" hPa ";
		//last_sub_waether_M<-"Coefficient de marée:"+tidal_coeff;
		
		write("weather after init : "+weather_info);
	}
}
/* ******************************************************************
 ******* fin global 			*******                                       ***
*********************************************************************/



/***************************************
 * ******* GRILLE DE PARCELLES ******* *
 ***************************************/

grid parcelle width: 52 height: 90 neighbors: 8 use_regular_agents: false use_individual_shapes: false 
{
	//type parcelle (0:mer, 1:nat, 2: agricole; 3: urb)
	int type_cell;
	
	/*************************** GESTION SUBMERSION ****************************************/
	float altitude;				// altitude d'apres le MNT
	bool is_sea <- false;		// cellule mer / terre 

	// hauteur d'eau sur la cellule
	float water_height <- 0.0 min:0.0;
	bool is_water<- false;
	float water_altitude <- 0.0; 
	float water_height_history <- 0.0; 
	
	// hauteur de vague sur la cellule
	float wave_height <- 0.0 min:0.0;
	float wave_altitude <-0.0;
	
	int nb_sub <- 0 ;		// valeur historique:submersion (nombre et hauteur max) 

	
	/******************** GESTION DIGUES **************************************************** */
	//digues
	int dyke<- 0;					// 0:pas de digue, 1:digue en rocher , 2:digue en BA et 3:dune
	float dyke_height <- 0.0;		// hauteur de la digue si digue
	
	//altitue digue = altitue+hauteur digue, sauf si elle est cassée et à ce moment altitude digue = altitude
	float dyke_altitude {
		float alti <-altitude+dyke_height;
		if bursted_dyke {alti<-altitude;}
		return alti;
	}
	int burst_coeff {
		if dyke=1 {return 3;}
		if dyke=2 {return 1;}
		if dyke=3 {return 1;}
	}
	
	float dyke_state;   //0 bursted, 100:état parfait
	bool bursted_dyke<-false;	//état de la digue (si state=0 alors busrted)
	
	//parcelle avec digue à coté
	list<parcelle> close_dyke;
	
			
	/************************************ GESTION CELLULES VOISINES ******************************* */
	// cellules voisines (Moore, 8)
	list<parcelle> neighbour_cells;
	list<parcelle> neighbour_cells_far;
	
	//parcelles de terre voisines (utilisé pour submersion)
	list<parcelle> earth_neighbour_cells;
	bool is_earth_neighbour_cells;
	
	//liste des parcelles de mer voisine de la parcelle
	list<parcelle> sea_neighbour_cells;
	bool close_sea <- false;
	float distance_to_sea;	// distance a la mer
	// parcelle de mer la plus proche et distance à la mer 
   	parcelle closestSea <- (sea_cells closest_to(self));
	
	// est-ce que la cellule a deja ete traitee dans la diffusion de l'eau
	bool already <- false;

	//ligne et colonne (pour les vagues)
	int line<-0;
	int row <-0;
	
	
	/********************************* INITIALISATION CELLULES ****************************************** */
	int init_parcel(matrix init_data, int pop) {
		altitude<-float(init_data[0,pop]);
		is_sea<-bool(init_data[1,pop]);
		water_height_history<-float(init_data[6,pop]);
			
		//pour les parcelles de mer (is_sea)
		if (is_sea) {
			water_height<-(-altitude);								//on démarre avec une altitude d'eau=0 et donc une hauteur d'eau=-altitude
			water_altitude <-0.0;
			is_water<-true;											//y'a de l'eau
		}
		
		dyke<-int(init_data[2,pop]);
		if dyke>0 {			//pour les parcelles avec digue
			dyke_height<-float(init_data[4,pop]) #m;
			is_sea<-false;
			water_height<-0.0;
			altitude<-max([0,altitude]);
			dyke_state<-90.0;
		} else {
			dyke_height<-0.0;
		}
			
		//définition des voisinages des parcelles
		neighbour_cells <- (self neighbors_at 1);
		neighbour_cells_far <- (self neighbors_at 2);
			
		// next column of data = next parcel
		return pop+1;
	}
	

	/****************************************************************************
	 * 			ACTIONS DE LA PARCELLE POUR LA SUBMERSION 	          *
	 ****************************************************************************/
	 
	//action deversement dans les parcelles voisines
	action flow {	
		//write("flow from cell "+name);
		ask earth_neighbour_cells	where (self.water_altitude> each.altitude)
		{
			// Test rupture de digue et altitude obstacle corrigée si besoin
			if (dyke>0) 
			{
				dyke_state<-max([0,dyke_state-myself.water_height*burst_coeff()]);
				if (dyke_state=0) 		{	bursted_dyke<-true;	}
			}
			
			if myself.water_altitude>dyke_altitude() 
			{
				water_altitude<-myself.water_altitude;
				water_height<- myself.water_altitude-altitude;
				is_water<- true;		//indique présence d'eau -> cette cellule pourra se déverser au tour prochain
			 }
		}
		already <- true;		//déjà déversé ce tour
	}//end action flow
		
	action wave_wind {
		//write("wave wind on cell "+name);
		//calcul hauteur de vague
		wave_altitude <- water_altitude+wave_height; 
		//on demande à la voisine de vent
		ask parcelle where (each.grid_x=(self.grid_x+wind_row) and each.grid_y=(self.grid_y+wind_line)) 
			{
				if (self.dyke>0) {
					dyke_state<-max([0,dyke_state-myself.wave_height*2*burst_coeff()]);
					if (dyke_state=0) 		{	bursted_dyke<-true;	}
				}//end if 1
			
				if (self.dyke_altitude()<myself.wave_altitude) 
				{
					if self.water_altitude>=myself.water_altitude 
					{
						self.wave_altitude<-myself.wave_altitude;
						self.is_water<-true;
					}
					else 
					{
						self.water_altitude<-self.altitude+myself.wave_altitude-self.dyke_altitude();
						self.is_water<-true;	
					}
				}//end if 2
		}// end ask
	}//end action
	
	//*************** DEFINITION DES CARTES DE VISUALISATION***************
	
	aspect sub_history	{
		// submersion history
		//if (!is_sea and water_height_history>0) {draw square(300 # m) color:rgb(130,200,255) border:#black;}

		//digue = triangle jaune (dune) ou gris(béton) ou marron (enrochement)
		if dyke>0 {
			if (dyke=1){draw triangle(100 # m) color:rgb(140,92,50) border:#black;}
			if (dyke=2){draw triangle(100 # m) color:rgb(180,180,180) border:#black; }
			if (dyke=3){draw triangle(100 # m) color:# yellow border:#black;}
		}
		
		// couleur fond de la cellule
		if (is_sea) {color<- #blue;}
		else if (!is_sea and is_water) {
			color<- rgb(116,208,241);
		}
		else if (!is_sea and !is_water) {
			color<-rgb(int(min([255,max([245 - 10 *altitude, 0])])), 
						int(min([255,max([245 - 15 *altitude, 0])])), 
						int(min([255,max([0,220 - 25 * altitude])]))
			);
		}
		else {
			color <- #orange;
		}
	}

}// end grid






/********************
 * *** SIMULATION ***
 ********************/
experiment Start type:gui
{
	//Definition de quelques parametres
	parameter "Submersion intensity: " var: submersion_intensity <- rnd(10)  min: 1 max: 10 category: "Meteo" ;
	parameter "Wind direction: " var: direction <- rnd(0,7) min: 0 max: 7 category: "Meteo";
	
	output
	{
	//Carte de la dernière submersion	
		display submersion name:"Submersion" ambient_light:100		
		{
			grid parcelle triangulation:false lines:# black;
			species parcelle aspect:sub_history;
		   	overlay position:{ 2, 2 } size:{ 300 #px, 150 #px } background:# black transparency:0.5 border:#black rounded:true 
		   	{
		   			//int last_sub_year<-2010;										--> sub_year
					//string last_sub_name<-"Xynthia";									--> storm_name[sub_nb]
					//string last_sub_waether_P<-"Pression atmosphèrique : 969 hPa"; 			--> atmosph_pressure
					//string last_sub_waether_M<-"Coefficient de marée : 102"; 				--> tidal_coeff
					//string last_sub_waether_V<-"Vent de Sud Ouest d'une force de 10 Bf";		--> weather_info
		   		
					//légende 
    		         	draw "Tempête : "+storm_name[sub_nb] at:{20#px,20#px } color:#white font:font("SansSerif", 14, #bold);
	    		        draw "Année : "+sub_year at:{10#px,40#px } color:#white font:font("SansSerif", 12, #bold);
	             	draw "Atmospheric pressure: "+atmosph_pressure at:{10#px,60#px } color:#white font:font("SansSerif", 12, #bold);
	             	draw "Tidal coeff: "+tidal_coeff at:{10#px,80#px } color:#white font:font("SansSerif", 12, #bold);
	             	draw wind_info at:{10#px,100#px } color:#white font:font("SansSerif", 12, #bold);
	             	draw wave_info at: {10#px,120#px } color:#white font:font("SansSerif", 12, #bold);
	             	draw water_info at: {10#px,140#px } color:#white font:font("SansSerif", 12, #bold);
			}	
		}// end display
	}//end output
}//end experiment 
//experiment Mode_test type:batch repeat:5 keep_seed:true until:(endcond=true){}

