

class timeData
{
 public:
    scalar t;
    scalar A;
    scalar z1;
	scalar z2;
	scalar alpha; 
    scalar U1;
    scalar U2;
    scalar vol;
    scalar viscP;
    scalar viscP_1;
	scalar Ca_1;
	scalar Ca_2;
        scalar p1;
        scalar p2;
    
	timeData()
	{
		t=0.0;
		A=0.0;
		z1=0.0;
		z2=0.0;
	    alpha=0.0; 
        U1=0.0;
        U2=0.0;
        vol=0.0;
        viscP=0.0;
        viscP_1=0.0;
		Ca_1= 0.0;
		Ca_2= 0.0;
	};
	    
};

class timeDataList 
{
 public:
    List<timeData> data;
 
    timeDataList(label size);
	
	timeDataList(const timeDataList& dataList, label size, label offset);
};

timeDataList::timeDataList(label size):
	data(size,timeData())
{
};

timeDataList::timeDataList(const timeDataList& inputDataList, label size, label offset):
	data(size,timeData())
{
	const List<timeData>& inputData = inputDataList.data;

	int size_group = (inputDataList.data.size()-offset)/size;
	double max  = 0.0;
	double min  = 0.0;
	double sum  = 0.0;
	double m    = 0.0;

    for (int i =0;i<size;i++)
	{
        int beg=i*size_group+offset;
		int end=(i+1)*size_group+offset;
		//t
		max  = inputData[beg].t;
		min  = inputData[beg].t;
		sum  = inputData[beg].t;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].t);
			min = std::min(min,inputData[j].t);
			sum += inputData[j].t;
		}
        m = sum/size_group;
        data[i].t=(sum+2*m-max-min)/size_group;

		//A
		max  = inputData[beg].A;
		min  = inputData[beg].A;
		sum  = inputData[beg].A;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].A);
			min = std::min(min,inputData[j].A);
			sum += inputData[j].A;
		}
        m = sum/size_group;
        data[i].A=(sum+2*m-max-min)/size_group;

		//z1
		max  = inputData[beg].z1;
		min  = inputData[beg].z1;
		sum  = inputData[beg].z1;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].z1);
			min = std::min(min,inputData[j].z1);
			sum += inputData[j].z1;
		}
        m = sum/size_group;
        data[i].z1=(sum+2*m-max-min)/size_group;

		//z2
		max  = inputData[beg].z2;
		min  = inputData[beg].z2;
		sum  = inputData[beg].z2;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].z2);
			min = std::min(min,inputData[j].z2);
			sum += inputData[j].z2;
		}
        m = sum/size_group;
        data[i].z2=(sum+2*m-max-min)/size_group;

		//alpha
		max  = inputData[beg].alpha;
		min  = inputData[beg].alpha;
		sum  = inputData[beg].alpha;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].alpha);
			min = std::min(min,inputData[j].alpha);
			sum += inputData[j].alpha;
		}
        m = sum/size_group;
        data[i].alpha=(sum+2*m-max-min)/size_group;

		//U1
		max  = inputData[beg].U1;
		min  = inputData[beg].U1;
		sum  = inputData[beg].U1;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].U1);
			min = std::min(min,inputData[j].U1);
			sum += inputData[j].U1;
		}
        m = sum/size_group;
        data[i].U1=(sum+2*m-max-min)/size_group;

		//U2
		max  = inputData[beg].U2;
		min  = inputData[beg].U2;
		sum  = inputData[beg].U2;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].U2);
			min = std::min(min,inputData[j].U2);
			sum += inputData[j].U2;
		}
        m = sum/size_group;
        data[i].U2=(sum+2*m-max-min)/size_group;

		//vol
		max  = inputData[beg].vol;
		min  = inputData[beg].vol;
		sum  = inputData[beg].vol;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].vol);
			min = std::min(min,inputData[j].vol);
			sum += inputData[j].vol;
		}
        m = sum/size_group;
        data[i].vol=(sum+2*m-max-min)/size_group;

		//vicsP
		max  = inputData[beg].viscP;
		min  = inputData[beg].viscP;
		sum  = inputData[beg].viscP;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].viscP);
			min = std::min(min,inputData[j].viscP);
			sum += inputData[j].viscP;
		}
        m = sum/size_group;
        data[i].viscP=(sum+2*m-max-min)/size_group;

		//viscP_1
		max  = inputData[beg].viscP_1;
		min  = inputData[beg].viscP_1;
		sum  = inputData[beg].viscP_1;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].viscP_1);
			min = std::min(min,inputData[j].viscP_1);
			sum += inputData[j].viscP_1;
		}
        m = sum/size_group;
        data[i].viscP_1=(sum+2*m-max-min)/size_group;

		//Ca_1
		max  = inputData[beg].Ca_1;
		min  = inputData[beg].Ca_1;
		sum  = inputData[beg].Ca_1;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].Ca_1);
			min = std::min(min,inputData[j].Ca_1);
			sum += inputData[j].Ca_1;
		}
        m = sum/size_group;
        data[i].Ca_1=(sum+2*m-max-min)/size_group;

		//Ca_2
		max  = inputData[beg].Ca_2;
		min  = inputData[beg].Ca_2;
		sum  = inputData[beg].Ca_2;
		for (int j=beg+1;j<end;j++)
		{
			max = std::max(max,inputData[j].Ca_2);
			min = std::min(min,inputData[j].Ca_2);
			sum += inputData[j].Ca_2;
		}
        m = sum/size_group;
        data[i].Ca_2=(sum+2*m-max-min)/size_group;
    }



};
