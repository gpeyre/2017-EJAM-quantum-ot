/*
 * (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
 * University, University College London. This code is released under the 
 * Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
 * license (http://creativecommons.org/licenses/by-nc-sa/4.0/).
 */

#define _CRT_SECURE_NO_WARNINGS
#include "App.hpp"
#include "base/Main.hpp"
#include "gpu/GLContext.hpp"
#include "3d/Mesh.hpp"
#include "io/File.hpp"
#include "io/StateDump.hpp"
#include "base/Random.hpp"

#include <stdio.h>
#include <conio.h>
#include <direct.h>
#include <stdlib.h>

using namespace FW;

//------------------------------------------------------------------------

App::App(void)
:   m_commonCtrl    (CommonControls::Feature_Default & ~CommonControls::Feature_RepaintOnF5),
    m_cameraCtrl    (&m_commonCtrl, CameraControls::Feature_Default | CameraControls::Feature_StereoControls),
    m_action        (Action_None),
	m_hasTextures	(false)
{
    m_commonCtrl.showFPS(false);
    m_commonCtrl.addStateObject(this);
    m_cameraCtrl.setKeepAligned(true);

	m_commonCtrl.addSlider(&m_intensity, 0.01f, 1000000.0f, true, FW_KEY_NONE, FW_KEY_NONE, "Intensity", 1.0f);
	m_commonCtrl.endSliderStack();
    m_commonCtrl.addSeparator();
	m_commonCtrl.beginSliderStack();
	m_commonCtrl.addSlider(&m_light_azimuth, 0.0f, 2*3.1415926f, false, FW_KEY_NONE, FW_KEY_NONE, "Light azimuth", 1.0f);
	m_commonCtrl.addSlider(&m_light_elevation, 0.0f, 0.5f*3.1415926f, false, FW_KEY_NONE, FW_KEY_NONE, "Light elevation", 1.0f);
	m_commonCtrl.addSlider(&m_light_distance, 0.01f, 1000.0f, true, FW_KEY_NONE, FW_KEY_NONE, "Light distance", 1.0f);
	m_commonCtrl.endSliderStack();
    m_window.addListener(&m_cameraCtrl);
    m_commonCtrl.addSeparator();

	m_commonCtrl.beginSliderStack();
	m_commonCtrl.addSlider(&m_obj_rot, 0.0f, 2*3.1415926f, false, FW_KEY_NONE, FW_KEY_NONE, "Plane rotation", 1.0f);
	m_commonCtrl.addSlider(&m_obj_scale, 1.0f, 1000.0f, true, FW_KEY_NONE, FW_KEY_NONE, "Plane scale", 1.0f);
	m_commonCtrl.endSliderStack();


	m_commonCtrl.beginSliderStack();
	m_commonCtrl.addSlider(&m_light_r, 0.0f, 1.0f, false, FW_KEY_NONE, FW_KEY_NONE, "Light color R", 1.0f);
	m_commonCtrl.addSlider(&m_light_g, 0.0f, 1.0f, false, FW_KEY_NONE, FW_KEY_NONE, "Light color G", 1.0f);
	m_commonCtrl.addSlider(&m_light_b, 0.0f, 1.0f, false, FW_KEY_NONE, FW_KEY_NONE, "Light color B", 1.0f);
	m_commonCtrl.endSliderStack();


	m_commonCtrl.addButton((S32*)&m_action, Action_LoadTex1,                FW_KEY_NONE,       "Load textures");

	m_window.setTitle("texrend");
    m_window.addListener(this);
    m_window.addListener(&m_commonCtrl);

	m_window.setSize(Vec2i(1920,1080));
	m_window.getGL();

	m_mtxCamera = m_cameraCtrl.getCameraToWorld();
	m_mtxObject.setIdentity();
	m_mtxEnv.setIdentity();
	m_controlMode = ControlMode_Camera;
	m_controlModePrev = ControlMode_None;

	// Default settings

	m_intensity = 10000.0f;
	m_light_azimuth = 3.0f;
	m_light_elevation = 3.14159f/4.0f;
	m_light_distance = 1.0f;

	/*
	// original flash position
	m_light_azimuth = 0.0f;
	m_light_elevation = 0.0f;
	m_light_distance = 100.0f;
	*/ 

	m_light_azimuth = 4.0f;
	m_light_elevation = 1.0f;
	m_light_distance = 55.0f;
	m_intensity = 1500.0f;

	m_obj_pos_x = m_obj_pos_y = 0;

	m_light_r = m_light_g = m_light_b = 1.0f;


	m_obj_scale = 100.0f;
	m_obj_rot = 0;

	m_hasTextures = true;
	m_use_fbo = false;
		

	initCamera();

	//loadMaterials(".\\temp\\asdf.pfm");

}



GLuint App::loadMaterials(String name)
{
	printf("Loading materials: %s\n", name.getDirName().getPtr());


	//// Maps
	m_texname1 = name.getDirName().append("\\map_diff.pfm");
	m_texname2 = name.getDirName().append("\\map_spec.pfm");
	m_texname3 = name.getDirName().append("\\map_spec_shape.pfm");
	m_texname4 = name.getDirName().append("\\map_normal.pfm");
	
	m_glTexture[1] = loadTexture(m_texname1);
	m_glTexture[2] = loadTexture(m_texname2);
	m_glTexture[3] = loadTexture(m_texname3);
	m_glTexture[7] = loadTexture(m_texname4);	
	

	m_tex_size = loadTextureSize(m_texname1);  // xxx

	//FILE* f = fopen("temp\\map_params.dat","rb");
	printf("%s\n", name.getDirName().append("\\map_params.dat").getPtr());
	FILE* f = fopen(name.getDirName().append("\\map_params.dat").getPtr(), "rb");
	fscanf(f, "%f", &m_alpha);
	printf("alpha %f\n", m_alpha);
	fclose(f);

	return 0;
}

GLuint App::loadTexture(String name)
{
	char tmp[256];	// xxx....
	char tmp2;
	int w, h;
	printf("Loading %s... ", name);
	FILE* f = fopen(name.getPtr(),"rb");
	fscanf(f, "%s\n%d %d\n-1.0%c", tmp, &w, &h, &tmp2);
	printf("(%i %s %i %i) ", ftell(f), tmp, w, h);
	Image img(Vec2i(w, h), ImageFormat::RGBA_Vec4f);
	fread(img.getMutablePtr(), 12, img.getSize().x * img.getSize().y, f);
	fclose(f);
	printf("%s\n%d %d\n-1.0\n", tmp, w, h);
	for (int i = img.getSize().x * img.getSize().y - 1; i >= 0; i--)
		((Vec4f*)img.getMutablePtr())[i] = Vec4f(((Vec3f*)img.getPtr())[i], 1.f); // vec3 -> vec4
	printf("Done\n");

	return img.createGLTexture();
}


Vec2f App::loadTextureSize(String name)
{
	Vec2f result;

	char tmp[256];
	char tmp2;
	int w, h;
	printf("Loading %s... ", name);
	FILE* f = fopen(name.getPtr(),"rb");
	fscanf(f, "%s\n%d %d\n-1.0%c", tmp, &w, &h, &tmp2);
	printf("SIZE: (%i %s %i %i) ", ftell(f), tmp, w, h);
	result.x = w;
	result.y = h;
	fclose(f);

	return result;
}



App::~App(void)
{
    // empty
}

//------------------------------------------------------------------------

String App::getTexImportFilter(void)
{
    return
		"pfm:Floating point";
}

String App::getPoseImportFilter(void)
{
    return
		"dat:Pose data";
}

bool App::handleEvent(const Window::Event& ev)
{
    if (ev.type == Window::EventType_Close)
    {
        m_window.showModalMessage("Exiting...");
        delete this;
        return true;
    }
#if 0
	// juggle matrices
	if (m_controlModePrev == ControlMode_Camera)
		m_mtxCamera = m_cameraCtrl.getWorldToCamera();
	else if (m_controlModePrev == ControlMode_Object)
		m_mtxObject = m_cameraCtrl.getWorldToCamera();
	else if (m_controlModePrev == ControlMode_Env)
	{
		m_mtxEnv = m_cameraCtrl.getWorldToCamera();
		Mat3f rot = m_mtxEnv.getXYZ();
		rot.col(1) = Vec3f(0.f, 1.f, 0.f);
		rot.col(0) = normalize(cross(rot.col(1), rot.col(2)));
		rot.col(2) = cross(rot.col(0), rot.col(1));
		m_mtxEnv.col(0) = Vec4f(rot.col(0), 0.f);
		m_mtxEnv.col(1) = Vec4f(rot.col(1), 0.f);
		m_mtxEnv.col(2) = Vec4f(rot.col(2), 0.f);
	}

	if (m_controlMode != m_controlModePrev)
	{
		if (m_controlMode == ControlMode_Camera)
			m_cameraCtrl.setWorldToCamera(m_mtxCamera);
		else if (m_controlMode == ControlMode_Object)
			m_cameraCtrl.setWorldToCamera(m_mtxObject);
		else if (m_controlMode == ControlMode_Env)
			m_cameraCtrl.setWorldToCamera(m_mtxEnv);
		m_controlModePrev = m_controlMode;
	}
#endif


    if (ev.type == Window::EventType_KeyDown)
    { 
		if (ev.key == FW_KEY_PAGE_UP)
			m_obj_scale *= 1.005f;
		if (ev.key == FW_KEY_PAGE_DOWN)
			m_obj_scale /= 1.005f;
		if (ev.key == FW_KEY_O)
			m_obj_rot = m_obj_rot+0.01f;
		if (ev.key == FW_KEY_U)
			m_obj_rot = m_obj_rot-0.01f;
		if (ev.key == FW_KEY_J)
			m_obj_pos_x -= 0.003*m_obj_scale;
		if (ev.key == FW_KEY_L)
			m_obj_pos_x += 0.003*m_obj_scale;
		if (ev.key == FW_KEY_K)
			m_obj_pos_y -= 0.003*m_obj_scale;
		if (ev.key == FW_KEY_I)
			m_obj_pos_y += 0.003*m_obj_scale;
		if (ev.key == FW_KEY_A)
			m_cameraCtrl.setUp(Vec3f(0,0,-100));
    }


    Action action = m_action;
    m_action = Action_None;
    String name;
    Mat4f mat;

    switch (action)
    {
    case Action_None:
        break;

    case Action_ResetCamera:
		initCamera();
        break;

    case Action_EncodeCameraSignature:
        m_window.setVisible(false);
        printf("\nCamera signature:\n");
        printf("%s\n", m_cameraCtrl.encodeSignature().getPtr());
        waitKey();
        break;

    case Action_DecodeCameraSignature:
        {
            m_window.setVisible(false);
            printf("\nEnter camera signature:\n");

            char buf[1024];
            if (scanf_s("%s", buf, FW_ARRAY_SIZE(buf)))
                m_cameraCtrl.decodeSignature(buf);
            else
                setError("Signature too long!");

            if (!hasError())
                printf("Done.\n\n");
            else
            {
                printf("Error: %s\n", getError().getPtr());
                clearError();
                waitKey();
            }
        }
		break;
	case Action_LoadTex1:
		{
			name = m_window.showFileLoadDialog("Load texture", getTexImportFilter());
			loadMaterials(name);

			printf("%s\n", name.getDirName().append("\\map_params2.dat").getPtr());
			FILE* f = fopen(name.getDirName().append("\\map_params2.dat").getPtr(), "rb");
			if (f)
			{
				float p1, p2, p3, p4, p5;
				fscanf(f, "%f %f %f %f %f", &p1, &p2, &p3, &p4, &p5);
				printf("params2 %f %f %f %f %f\n", p1, p2, p3, p4, p5);
				fclose(f);
			}
		}		
        break;

    default:
        FW_ASSERT(false);
        break;
    }

    m_window.setVisible(true);

    if (ev.type == Window::EventType_Paint)
        renderFrame(m_window.getGL());
    m_window.repaint();
    return false;
}


//------------------------------------------------------------------------

void App::readState(StateDump& d)
{
	// XXX not up to date
    d.pushOwner("App");
	d.get(m_mtxCamera,         "m_mtxCamera");
	d.get(m_mtxObject,         "m_mtxObject");
	d.get(m_mtxEnv,            "m_mtxEnv");
	d.get((S32&)m_controlMode, "m_controlMode");
	m_controlModePrev = ControlMode_None;

	d.get(m_intensity,			"m_intensity");
	d.get(m_light_azimuth,		"m_light_azimuth");
	d.get(m_light_elevation,	"m_light_elevation");
	d.get(m_light_distance,		"m_light_distance");

	d.get(m_texname1,			"m_texname1");
	d.get(m_texname2,			"m_texname2");
	d.get(m_texname3,			"m_texname3");
	d.get(m_texname4,			"m_texname4");

	glDeleteTextures(1,&m_glTexture[1]);
	m_glTexture[1] = loadTexture(m_texname1);
	glDeleteTextures(1,&m_glTexture[2]);
	m_glTexture[2] = loadTexture(m_texname2);
	glDeleteTextures(1,&m_glTexture[3]);
	m_glTexture[3] = loadTexture(m_texname3);
	glDeleteTextures(1,&m_glTexture[7]);
	m_glTexture[7] = loadTexture(m_texname4);

    d.popOwner();
}

//------------------------------------------------------------------------

void App::writeState(StateDump& d) const
{
    d.pushOwner("App");
	d.set(m_mtxCamera,        "m_mtxCamera");
	d.set(m_mtxObject,        "m_mtxObject");
	d.set(m_mtxEnv,           "m_mtxEnv");
	d.set((S32)m_controlMode, "m_controlMode");

	d.set(m_intensity,			"m_intensity");
	d.set(m_light_azimuth,		"m_light_azimuth");
	d.set(m_light_elevation,		"m_light_elevation");
	d.set(m_light_distance,		"m_light_distance");
	d.set(m_texname1,           "m_texname1");
	d.set(m_texname2,           "m_texname2");
	d.set(m_texname3,           "m_texname3");
	d.set(m_texname4,           "m_texname4");
    d.popOwner();
}

//------------------------------------------------------------------------

void App::waitKey(void)
{
    printf("Press any key to continue . . . ");
    _getch();
    printf("\n\n");
}



//------------------------------------------------------------------------

void App::renderFrame(GLContext* gl)
{	
	if (!m_hasTextures)
	{
        gl->drawModalMessage("No maps loaded.");
		return;
	}


    Mat4f worldToCamera = m_cameraCtrl.getWorldToCamera();


	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

	Mat4f projection = gl->xformFitToView(-Vec2f(-1.0f, -1.0f), -Vec2f(2.0f, 2.0f)) * m_cameraCtrl.getCameraToClip();

    renderQuad(gl, projection);	

}

//------------------------------------------------------------------------

void App::renderBackground(GLContext* gl, const Mat4f& projection)
{
    FW_ASSERT(gl);
    const char* progId = "backgroundshader";
    GLContext::Program* prog = gl->getProgram(progId);

    if (!prog)
    {
        prog = new GLContext::Program(
            "#version 120\n"
            FW_GL_SHADER_SOURCE(
                attribute vec2 xyIn;
                attribute vec3 dirIn;
                varying vec3 dir;
				

                void main()
                {
                    gl_Position = vec4(xyIn, 0.0, 1.0);
                    dir = dirIn;

					// XXX overwrite for now, 10.1.2013:
					dir.xy = (xyIn*0.5+0.5)*vec2(1.0,-1.0);
                }
            ),
            "#version 120\n"
            FW_GL_SHADER_SOURCE(
                uniform sampler2D tex0;
                uniform sampler2D tex1;
                uniform sampler2D tex2;
                uniform sampler2D tex3;
                uniform sampler2D tex4;
				
                varying vec3 dir;

                void main()
                {
					vec3 ndir = normalize(dir);

					float PI = 3.14159265358979323846264;

					float u = (1+atan(ndir.x,-ndir.z)/PI)*0.5;
					float v = (acos(ndir.y)/PI);

					gl_FragColor = vec4(0);
                }
            ));
    }

	glPushAttrib(GL_ALL_ATTRIB_BITS);

	// Render the quad.

    gl->setProgram(progId, prog);
	prog->use();
    gl->setUniform(prog->getUniformLoc("tex0"), 0);
    gl->setUniform(prog->getUniformLoc("tex1"), 1);
    gl->setUniform(prog->getUniformLoc("tex2"), 2);
    gl->setUniform(prog->getUniformLoc("tex3"), 3);
    gl->setUniform(prog->getUniformLoc("tex4"), 4);
    gl->setUniform(prog->getUniformLoc("tex5"), 5);

	glEnableVertexAttribArray(prog->getAttribLoc("xyIn"));
	glEnableVertexAttribArray(prog->getAttribLoc("dirIn"));	


	for (int i=0; i < NTEXTURES; i++)
	{
		glActiveTexture(GL_TEXTURE0 + i);
		glBindTexture(GL_TEXTURE_2D, m_glTexture[i]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, FW_S32_MAX);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

	}

	Vec2f xy[4] = {
		Vec2f(-1.f, -1.f),
		Vec2f(1.f, -1.f),
		Vec2f(1.f, 1.f),
		Vec2f(-1.f, 1.f)
	};

	Vec3f dir[4];
	for (int i=0; i < 4; i++)
	{
		Vec4f p(xy[i], 0.f, 1.f);
		p = invert(projection) * p;
		dir[i] = invert(m_mtxCamera.getXYZ()) * p.toCartesian();
		dir[i] = m_mtxEnv * dir[i];
	}
	
	glVertexAttribPointer(prog->getAttribLoc("xyIn"),  2, GL_FLOAT, false, 0, xy);
	glVertexAttribPointer(prog->getAttribLoc("dirIn"), 3, GL_FLOAT, false, 0, dir);
	glDepthMask(GL_FALSE);
	glDrawArrays(GL_QUADS, 0, 4);

	gl->checkErrors();

	glPopAttrib();
}

//------------------------------------------------------------------------

void App::renderQuad(GLContext* gl, const Mat4f& projection)
{
    FW_ASSERT(gl);
    const char* progId = "quadshader";
    GLContext::Program* prog = gl->getProgram(progId);

    if (!prog)
    {
        prog = new GLContext::Program(
            "#version 120\n"
            FW_GL_SHADER_SOURCE(
                uniform mat4 worldToClip;
                uniform mat4 worldToCamera;
                attribute vec3 posIn;
                attribute vec2 uvIn;
				varying vec3 pos;
                varying vec2 uv;
				uniform float obj_rot;

				uniform vec2 obj_pos;
				


				uniform mat3 gt_T;
				uniform vec3 gt_c;

                void main()
                {
					vec4 pos4 = vec4(posIn, 1.0);

					float or = obj_rot;
					mat2 OR = mat2(cos(or), -sin(or),
							sin(or), cos(or));
					vec4 pos4_2 = pos4;
					pos4_2.xy = OR*pos4_2.xy;
					pos4_2.xy += obj_pos;
					
                    gl_Position = worldToClip * pos4_2;

					pos = pos4.xyz;
                    uv  = uvIn;
                }
            ),
            "#version 120\n"
            FW_GL_SHADER_SOURCE(
                uniform sampler2D tex0;
                uniform sampler2D tex1;
                uniform sampler2D tex2;
                uniform sampler2D tex3;
				uniform sampler2D tex4;
				uniform sampler2D tex5;
				uniform sampler2D tex6;
				uniform vec3 cameraPos;
				uniform mat3 dirToEnv;
				uniform vec3 normal;
                varying vec3 pos;
                varying vec2 uv;

				uniform float alpha;

				uniform float intensity;
				uniform float light_azimuth;
				uniform float light_elevation;
				uniform float light_distance;
				uniform vec3 light_color;

				uniform float obj_rot;
				uniform vec2 obj_pos;

                void main()
                {

					vec3 alb_d = texture2D(tex1, uv).rgb;
					vec3 alb_s = texture2D(tex2, uv).rgb;
					vec3 specv = texture2D(tex3, uv).rgb;

					float or = -obj_rot;
					mat2 OR = mat2(cos(or), -sin(or),
							sin(or), cos(or));

					vec3 pos2 = pos;
					vec3 cameraPos2 = cameraPos;
					cameraPos2.xy = OR*(cameraPos2.xy-obj_pos);
					vec3 E = normalize(cameraPos2-pos2);

					vec3 L;

					L.x = sin(light_elevation)*cos(light_azimuth);
					L.y = sin(light_elevation)*sin(light_azimuth);
					L.z = cos(light_elevation);
					L = L * light_distance;
					L.xy = OR * (L.xy-obj_pos);

					L = L - pos2;


					float D2 = dot(L,L);

					L = L / sqrt(D2);


					float h = light_distance*cos(light_elevation);

					vec3 H = normalize(L+E);


					vec3 N = texture2D(tex4, uv).rgb;
					N = N / length(N);


					mat3 R = mat3(
						0, 0, N.x,
						0, 0, N.y,
						-N.x, -N.y, 0);

					// Halfway vector in normal-oriented coordinates (so normal is [0,0,1])
					vec3 Hn = H + R*H + 1.0/(N.z+1.0) * (R*(R*H));

					vec3 Hnp = Hn / vec3(Hn.z);
					
					mat2 M = mat2(
						specv.x, specv.z,
						specv.z, specv.y);

					vec3 HnpW = vec3(M * Hnp.xy, 1.0);

					float spec = exp(-pow(dot(HnpW.xy,Hnp.xy), alpha*0.5));	

					float cosine = max(.0, dot(N,L));

					float F0 = 0.04;

					float fres = F0 + (1-F0)*pow(1.0-max(0,(dot(H,E))), 5.0);	// Schlick
					spec = spec * fres / F0;
					
					spec = spec / vec3(dot(H, L));		// From Brady et al. model A

					gl_FragColor.rgb = (vec3(spec) * alb_s + alb_d)*vec3(cosine)/D2 * intensity * light_color;

					gl_FragColor.rgb = sqrt(gl_FragColor.rgb);	// rough gamma
					gl_FragColor.a = 1.0;

                }
            ));
    }

	glPushAttrib(GL_ALL_ATTRIB_BITS);

	// Render the quad.

    gl->setProgram(progId, prog);
	prog->use();

	m_mtxCamera = m_cameraCtrl.getWorldToCamera();


    gl->setUniform(prog->getUniformLoc("worldToClip"), projection * m_mtxCamera);
    gl->setUniform(prog->getUniformLoc("worldToCamera"), m_mtxCamera);
    gl->setUniform(prog->getUniformLoc("dirToEnv"), m_mtxEnv.getXYZ());

    gl->setUniform(prog->getUniformLoc("cameraPos"), Vec4f(invert(m_mtxCamera).getCol(3)).getXYZ());
    gl->setUniform(prog->getUniformLoc("normal"), m_mtxObject.getXYZ() * Vec3f(0.f, 1.f, 0.f));
    gl->setUniform(prog->getUniformLoc("tex0"), 0);
    gl->setUniform(prog->getUniformLoc("tex1"), 1);
    gl->setUniform(prog->getUniformLoc("tex2"), 2);
    gl->setUniform(prog->getUniformLoc("tex3"), 3);
	gl->setUniform(prog->getUniformLoc("tex4"), 7);
	gl->setUniform(prog->getUniformLoc("tex5"), 5);
	gl->setUniform(prog->getUniformLoc("tex6"), 6);
	glEnableVertexAttribArray(prog->getAttribLoc("posIn"));

	glEnableVertexAttribArray(prog->getAttribLoc("uvIn"));


	gl->setUniform(prog->getUniformLoc("alpha"), m_alpha);
	
	gl->setUniform(prog->getUniformLoc("intensity"), m_intensity);
	gl->setUniform(prog->getUniformLoc("light_azimuth"), m_light_azimuth);
	gl->setUniform(prog->getUniformLoc("light_elevation"), m_light_elevation);
	gl->setUniform(prog->getUniformLoc("light_distance"), m_light_distance);

	gl->setUniform(prog->getUniformLoc("light_color"), Vec3f(m_light_r,m_light_g,m_light_b));

	gl->setUniform(prog->getUniformLoc("obj_rot"), m_obj_rot);

	gl->setUniform(prog->getUniformLoc("divide"), (GetAsyncKeyState(VK_SPACE) & 0x8000));

	gl->setUniform(prog->getUniformLoc("win_size"), (Vec2f)m_window.getSize());

	gl->setUniform(prog->getUniformLoc("obj_pos"), Vec2f(m_obj_pos_x, m_obj_pos_y));
	

	for (int i=0; i < NTEXTURES; i++)
	{
		glActiveTexture(GL_TEXTURE0 + i);
		glBindTexture(GL_TEXTURE_2D, m_glTexture[i]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, FW_S32_MAX);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

		
		if (i == 4) {
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		}
	}

	Mat3f R;
	R.setZero();
	R(0,0) = cos(m_obj_rot);
	R(0,1) = -sin(m_obj_rot);
	R(1,0) = sin(m_obj_rot);
	R(1,1) = cos(m_obj_rot);
	R(2,2) = 1.0f;

	float aspect = (float)m_tex_size.y / m_tex_size.x;


	const float planesize = m_obj_scale*0.5f;
	Vec3f pos[4] = {
		Vec3f(-planesize, -planesize*aspect, 0.f),
		Vec3f(planesize, -planesize*aspect, 0.f),
		Vec3f(planesize, planesize*aspect, 0.f),
		Vec3f(-planesize, planesize*aspect, 0.f)
	};

	Vec2f uv[4] = {
		Vec2f(0.f, 0.f),
		Vec2f(1.f, 0.f),
		Vec2f(1.f, 1.f),
		Vec2f(0.f, 1.f)
	};

	glVertexAttribPointer(prog->getAttribLoc("posIn"), 3, GL_FLOAT, false, 0, pos);
	glVertexAttribPointer(prog->getAttribLoc("uvIn"), 2, GL_FLOAT, false, 0, uv);
	glDrawArrays(GL_QUADS, 0, 4);

	gl->checkErrors();

	glPopAttrib();
}

//------------------------------------------------------------------------

void App::initCamera(void)
{
    m_cameraCtrl.initDefaults();
	m_cameraCtrl.setKeepAligned(true);
	m_cameraCtrl.setNear(.1f);
	m_cameraCtrl.setFar(1000.f);

	/*
	// original head-on input photo position
	m_cameraCtrl.setPosition(Vec3f(0.1,0.1,100));
	m_cameraCtrl.setForward(Vec3f(0,0,-100));
	m_cameraCtrl.setUp(Vec3f(0,-100,0));
	m_cameraCtrl.setFOV(55);
	*/

	Vec3f lookat = Vec3f(0,10,-0.0f);
	Vec3f pos = Vec3f(5.0f,55.0f,60.0f);
	Vec3f forward = lookat-pos;
	m_cameraCtrl.setPosition(pos);
	m_cameraCtrl.setForward(forward);
	m_cameraCtrl.setUp(Vec3f(0,0,-1));
	m_cameraCtrl.setFOV(75);


    m_commonCtrl.message("Camera reset");
}

//------------------------------------------------------------------------

void FW::init(void)
{
    new App;
}

//------------------------------------------------------------------------

