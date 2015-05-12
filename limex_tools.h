#ifndef LIMEX_TOOLS_H__
#define LIMEX_TOOLS_H__

#define PROFILE_LIMEX
#ifdef PROFILE_LIMEX
	#define LIMEX_PROFILE_FUNC()		PROFILE_FUNC_GROUP("lime")
	#define LIMEX_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "limex")
	#define LIMEX_PROFILE_END()		PROFILE_END()
#else
	#define LIMEX_PROFILE_FUNC()
	#define LIMEX_PROFILE_BEGIN(name)
	#define LIMEX_PROFILE_END()
#endif

endif
