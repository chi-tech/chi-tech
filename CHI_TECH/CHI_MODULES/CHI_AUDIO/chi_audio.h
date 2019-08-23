#ifndef CHI_AUDIO_H
#define CHI_AUDIO_H

#include "../../CHI_VECTOR/chi_vector.h"
//#include "CHI_SOUND/chi_sound.h"
#include "chi_audio_incdef.h"

#include "../../../CHI_RESOURCES/Dependencies/stb/stb_vorbis.h"
#include "../../../CHI_RESOURCES/Dependencies/tiny/tinysound.h"

// template class CHI_VECTOR<CHI_SOUND>;

//#############################################################################
//Class definition
/***/
/** Object for handling Audio.
\author GMO */
class CHI_AUDIO
{
    public:
    //==========================================-========== Common attributes
    //	CHI_VECTOR<CHI_SOUND>                       soundStack;					///< Sound Stack
    //	CHI_VECTOR<CHI_AUDIO>						audioStack;					///< Audio Stack
    //
    int 										frequency;					///< Frequency of audio
    int 										latency;					///< Latency in Hz
    int 										bufSeconds;					///< Time buffer will hold
    int 										playingPool;
    int 										numPlayingPool;				///< Array for playing sounds

    char 										name[1000];					///< Filename

    tsLoadedSound 								vorbisLoaded;				///< Loaded OGG File
    tsPlaySoundDef								vorbisSound;				///< Playing OGG File
    tsContext*									context;					///< Initialize

    public:
    //======================-============================== Methods
    // 00 ConstrDestr
    CHI_AUDIO();
    ~CHI_AUDIO();

    // 01 Initialize
   	bool					LoadVorbis(char* fileName);
   	bool 					CleanAudio();

   	//02 Controls
   	void 					PlayAudio(int loopAmount);
   	void 					StopAudio();
   	void 					VolumeAudio(float volume);
};

#endif
