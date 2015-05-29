#include <vector>

class track
{
	public:
   float Momentum;  
   float pt;        
   float ptError;   
   float phi;			
   float eta;			
	float vx,vy,vz;	
	float dEdx;
	int charge;
	bool isKaon;
	bool vflag;
};

class compositeParticle
{
   public:
   float m;
   float eta;
   float phi;
   float pt; 
   float pt1;
   float pt2;
   float phi1;
   float phi2;
   float eta1;
   float eta2;
	float dEdx1;
	float dEdx2;
   float ptError1;
   float ptError2;
	bool vflag;
};


class EventData 
{

	public:
	int EventID;
   std::vector<track> tracks;
   std::vector<compositeParticle> compositeParticles;	
	float zVtx; 
	int multiplicityClass;

	void AddTrack(const track& p) { tracks.push_back(p); }
	int GetnTracks () { return tracks.size(); }

	void AddCompositeParticle(const compositeParticle& p) { compositeParticles.push_back(p); }
	int GetncompositeParticles () {return compositeParticles.size(); };

	void SetzVtx(float zVtx_) { zVtx = zVtx_; }
	void SetMultiplicityClass(int mclass_) { multiplicityClass  = mclass_; }

	void Clear() { tracks.clear(); compositeParticles.clear();}
};
