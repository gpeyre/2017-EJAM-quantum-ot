/*
 * (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
 * University, University College London. This code is released under the 
 * Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
 * license (http://creativecommons.org/licenses/by-nc-sa/4.0/).
 */

#pragma once
#include "gui/Window.hpp"
#include "gui/CommonControls.hpp"
#include "3d/CameraControls.hpp"
#include "gpu/Buffer.hpp"

namespace FW
{
//------------------------------------------------------------------------

class App : public Window::Listener, public CommonControls::StateObject
{
private:
    enum Action
    {
        Action_None,

        Action_ResetCamera,
        Action_EncodeCameraSignature,
        Action_DecodeCameraSignature,
		Action_LoadTex1,
		Action_LoadTex2,
		Action_LoadTex3,
		Action_LoadTex4,
    };

    enum ControlMode
    {
        ControlMode_None = 0,
		ControlMode_Camera,
		ControlMode_Object,
		ControlMode_Env
    };




public:
                    App             (void);
    virtual         ~App            (void);

    virtual bool    handleEvent     (const Window::Event& ev);
    virtual void    readState       (StateDump& d);
    virtual void    writeState      (StateDump& d) const;

private:
    void            waitKey         (void);
    void            renderFrame     (GLContext* gl);
	void			renderBackground(GLContext* gl, const Mat4f& projection);
    void            renderQuad      (GLContext* gl, const Mat4f& projection);

	Vec2f			loadTextureSize	(String name);
	GLuint			loadTexture		(String name);
	GLuint			loadTexturePNG		(String name, int mapping = 0);
	GLuint			loadMaterials	(String name);

	String getTexImportFilter(void);
	String getPoseImportFilter(void);

private:
                    App             (const App&); // forbidden
    App&            operator=       (const App&); // forbidden

	void			initCamera		(void);

private:
    Window          m_window;
    CommonControls  m_commonCtrl;
    CameraControls  m_cameraCtrl;

	Mat4f			m_mtxCamera;
	Mat4f			m_mtxObject;
	Mat4f			m_mtxEnv;
	ControlMode		m_controlMode;
	ControlMode		m_controlModePrev;

    Action          m_action;
    String          m_meshFileName;

	const static int NTEXTURES = 8;
	bool			m_hasTextures;
	U32				m_glTexture[NTEXTURES];

	F32				m_intensity;
	F32				m_light_azimuth;
	F32				m_light_elevation;
	F32				m_light_distance;
	F32				m_light_r;
	F32				m_light_g;
	F32				m_light_b;

	F32				m_obj_rot;
	F32				m_obj_scale;
	F32				m_obj_pos_x;
	F32				m_obj_pos_y;

	F32				m_alpha;

	String			m_texname1;
	String			m_texname2;
	String			m_texname3;
	String			m_texname4;
	String			m_texname5;
	String			m_texname6;

	Vec2f			m_tex_size;

	bool				m_use_fbo;
    GLuint              m_glDrawColorRenderbuffer;
    GLuint              m_glDrawDepthRenderbuffer;
	GLuint				m_glResolveColorRenderbuffer;
    GLuint              m_glDrawFramebuffer;
	GLuint				m_glResolveFramebuffer;
	Vec2i				m_resolution;
};



//------------------------------------------------------------------------
}
