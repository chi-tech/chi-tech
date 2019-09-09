#include "chi_thermoalpha.h"
#include "hydro_boundary/hydro_boundary.h"
#include "hydro_volume/hydro_volume.h"

//######################################################### Initialize
/** Initializes a thermal-hydraulic system.
 *
 *
 * \return Returns true when success.
\author Jan*/
bool CHI_THERMOSYSTEM::Initialize()
{
	//=========================================== Assigning fluid properties
	switch (fluidOption)
	{
		default:
		{
			CHI_FLOWPROPERTIES_WATER* waterProps = new CHI_FLOWPROPERTIES_WATER;
			flProps = (CHI_FLOWPROPERTIES*)waterProps;
			break;
		}
	}
	
	//=========================================== Open output file
	FILE* outpFile = fopen("ZZOUTPUT.txt","w");
	if (outpFile==NULL)
	{
		printf("ERROR: Could not open %s.\n","ZZOUTPUT.txt");
		return 0;
	}
	else
	{
		fprintf(outpFile,"################################# THERMOALPHA #################################\n\n");
	}
	
	
	//=========================================== Create traversal record
	CHI_VECTOR<bool> traverseRecord;
	for (int k=0;k<hydroComponentStack.itemCount;k++)
	{
		bool* newBool = new bool; *newBool = false;
		traverseRecord.PushItem(newBool);
	}
	
	//=========================================== Find first BC
	HYDRO_COMPONENT* ltComponent=NULL;
	for (int k=0;k<hydroComponentStack.itemCount;k++)
	{
		HYDRO_COMPONENT* curComp = hydroComponentStack.GetItem(k);
		if (curComp!=NULL)
		{

			if (curComp->componentType==2)
			{
				ltComponent = curComp;
				nodalization.PushItem(ltComponent);
				break;
			}
		}
	}
	
	
	//=========================================== Traversing
	HYDRO_COMPONENT* rtComponent;
	bool endTraverse=false;
	bool initialTravers = true;
	fprintf(outpFile,"\n=========================================================== Connectivity\n");
	while (!endTraverse)
	{
		
		int ltIndex = -1;
		int rtIndex = -1;
		int jcIndex = -1;
		
		//=================== Find junction index
		ltIndex = ltComponent->index;
		if (ltComponent->componentType==2)
		{
			HYDRO_BOUNDARY* bound = (HYDRO_BOUNDARY*)ltComponent;
			jcIndex = bound->junctionIndex;
		}
		else
		{
			HYDRO_VOLUME* volum = (HYDRO_VOLUME*)ltComponent;
			jcIndex = volum->rigtJunctionIndex;
		}

		
		//=================== Process junction
		HYDRO_SJUNCTION* junction = hydroSJunctionStack.GetItem(jcIndex);
		if (junction==NULL){break;}
		rtComponent = junction->rigtHComponent;
		if (rtComponent==NULL){break;}
		
		//=================== Find destination index
		rtIndex = rtComponent->index;
		nodalization.PushItem(rtComponent);
		
		
		fprintf(outpFile,"Component %03d, connected to Component %03d, via Junction %03d\n",ltIndex,rtIndex,jcIndex);
		
		if (rtComponent->componentType==2)
		{break;}
		ltComponent=rtComponent;
		
	}
	
	//=========================================== Initializing the nodalization
	fprintf(outpFile,"\n=========================================================== Initialized Parameters\n");
	fprintf(outpFile,"Node# Index Pressure    Tf       Tg       alpha_g Psat      rho_f  \n");
	for (int k=0;k<nodalization.itemCount;k++)
	{
		HYDRO_COMPONENT* curComp = nodalization.GetItem(k);
		if (curComp->componentType==1)
		{
			HYDRO_VOLUME* curVolume = (HYDRO_VOLUME*)curComp;
			double Psat = flProps->SaturationPressureT(curVolume->T_f);
			
			if (curVolume->P>Psat)
			{
				curVolume->rho_f = flProps->LiqDensityPT(curVolume->P,curVolume->T_f);
			}
			fprintf(outpFile,"%3d   %03d   %9.2f  %8.3f %8.3f %7.4f %9.2f %8.3f\n",k,
										curVolume->index,
										curVolume->P,
										curVolume->T_f,
										curVolume->T_g,
										curVolume->alpha_g,Psat,
										curVolume->rho_f);
		}
		if (curComp->componentType==2)
		{
			HYDRO_BOUNDARY* curVolume = (HYDRO_BOUNDARY*)curComp;
			double Psat = flProps->SaturationPressureT(curVolume->T_f);
			
			if (curVolume->P>Psat)
			{
				curVolume->rho_f = flProps->LiqDensityPT(curVolume->P,curVolume->T_f);
			}
			fprintf(outpFile,"%3d   %03d   %9.2f  %8.3f %8.3f %7.4f %9.2f %8.3f\n",k,
			        curVolume->index,
			        curVolume->P,
			        curVolume->T_f,
			        curVolume->T_g,
			        curVolume->alpha_g,Psat,
			        curVolume->rho_f);
		}
	}
	
	fclose(outpFile);
	
	//=========================================== Clearing traversal record
	traverseRecord.ClearVector();
	
	return false;
}