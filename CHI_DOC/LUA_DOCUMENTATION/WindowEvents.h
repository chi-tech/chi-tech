/**\defgroup LuaEvents Window Events
 * \ingroup LuaGeneralUtilities*/

/**Event generated when mouse moves
 * \ingroup LuaEvents*/
struct WM_MOUSEMOVE
{
    bool occured; ///< Flag indicating whether event occured
    int iPar0; ///< Mouse x-position
    int iPar1; ///< Mouse y-position
    int iPar2; ///< Mouse old x-position
    int iPar3; ///< Mouse old y-position
    int iPar4; ///< Scene/Window number
    int iPar5; ///< Object index under mouse pointer

    char sPar0; ///< Object name under mouse pointer
};

/**Event generated when left mouse button goes down
 * \ingroup LuaEvents*/
struct WM_LBUTTONDOWN
{
    bool occured; ///< Flag indicating whether event occured
    int iPar0; ///< Mouse x-position
    int iPar1; ///< Mouse y-position

    int iPar4; ///< Scene/Window number
    int iPar5; ///< Object index under mouse pointer

    char sPar0; ///< Object name under mouse pointer
};


/**Event generated when middle mouse button goes down
 * \ingroup LuaEvents*/
struct WM_MBUTTONDOWN
{
    bool occured; ///< Flag indicating whether event occured
    int iPar0; ///< Mouse x-position
    int iPar1; ///< Mouse y-position

    int iPar4; ///< Scene/Window number
    int iPar5; ///< Object index under mouse pointer

    char sPar0; ///< Object name under mouse pointer
};

/**Event generated when right mouse button goes down
 * \ingroup LuaEvents*/
struct WM_RBUTTONDOWN
{
    bool occured; ///< Flag indicating whether event occured
    int iPar0; ///< Mouse x-position
    int iPar1; ///< Mouse y-position

    int iPar4; ///< Scene/Window number
    int iPar5; ///< Object index under mouse pointer

    char sPar0; ///< Object name under mouse pointer
};

/**Event generated when left mouse button goes up
 * \ingroup LuaEvents*/
struct WM_LBUTTONUP
{
    bool occured; ///< Flag indicating whether event occured
    int iPar0; ///< Mouse x-position
    int iPar1; ///< Mouse y-position

    int iPar4; ///< Scene/Window number
    int iPar5; ///< Object index under mouse pointer

    char sPar0; ///< Object name under mouse pointer
};


/**Event generated when middle mouse button goes up
 * \ingroup LuaEvents*/
struct WM_MBUTTONUP
{
    bool occured; ///< Flag indicating whether event occured
    int iPar0; ///< Mouse x-position
    int iPar1; ///< Mouse y-position

    int iPar4; ///< Scene/Window number
    int iPar5; ///< Object index under mouse pointer

    char sPar0; ///< Object name under mouse pointer
};

/**Event generated when right mouse button goes up
 * \ingroup LuaEvents*/
struct WM_RBUTTONUP
{
    bool occured; ///< Flag indicating whether event occured
    int iPar0; ///< Mouse x-position
    int iPar1; ///< Mouse y-position

    int iPar4; ///< Scene/Window number
    int iPar5; ///< Object index under mouse pointer

    char sPar0; ///< Object name under mouse pointer
};

/**Event generated when mouse wheel has been rotated
 * \ingroup LuaEvents*/
struct WM_MOUSEWHEEL
{
    bool occured;
    int iPar0; ///< Wheel rotation (+-120)

    int iPar4; ///< Scene/Window number
    int iPar5; ///< Object index under mouse pointer

    char sPar0; ///< Object name under mouse pointer
};

/**Event generated when window is resized
 * \ingroup LuaEvents*/
struct WM_SIZE
{
    bool occured; ///< Flag indicating whether event occured

    int iPar4; ///< Scene/Window number

};

/**Event generated when keyboard key goes down
 * \ingroup LuaEvents*/
struct WM_KEYDN
{
    bool occured; ///< Flag indicating whether event occured
    int iPar0; ///< ASCII character code

    int iPar4; ///< Scene/Window number

};

/**Event generated when keyboard key goes up
 * \ingroup LuaEvents*/
struct WM_KEYUP
{
    bool occured; ///< Flag indicating whether event occured
    int iPar0; ///< ASCII character code

    int iPar4; ///< Scene/Window number

};

