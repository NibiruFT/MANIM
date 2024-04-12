from manim import *
from manim_voiceover import VoiceoverScene
from manim_voiceover.services.recorder import RecorderService

class Phasors(VoiceoverScene,MovingCameraScene):
    def construct(self):
        template = TexTemplate()
        template.add_to_preamble(r"\usepackage{tikz}\usepackage{circuitikz}\usepackage{gensymb}")
        Title = Tex(r'9.4 Phasor relationships for circuit elements', color=BLUE, font_size=50)
        self.set_speech_service(RecorderService(silence_threshold=-40.0))
        
        Resistor_Title = Tex(r'Resistor').move_to(RIGHT*4).scale(0.8) 
        Inductor_Title = Tex(r'Inductor').move_to(RIGHT*4).scale(0.8) 
        Capacitor_Title = Tex(r'Capacitor').move_to(RIGHT*4).scale(0.8) 
        frameboxRES = SurroundingRectangle(Resistor_Title, buff = .1)
        frameboxIND = SurroundingRectangle(Inductor_Title, buff = .1)
        frameboxCAP = SurroundingRectangle(Capacitor_Title, buff = .1)
    # Sine and Cosine and shifted
        sin_func = FunctionGraph(lambda t: np.sin(t) + 0.5 ,color=BLUE)
        sin_func.add_updater(lambda mob,dt: mob.become(FunctionGraph(lambda t: np.sin(t-dt) + 0.5, x_range =[-2*PI, 2*PI+self.renderer.time+dt] ,color=BLUE)))
        sin_func.add_updater(lambda mob,dt: mob.become(mob.shift(LEFT*+(self.renderer.time+dt-1))))
        
        cos_func = FunctionGraph(lambda t: np.cos(t) + 0.5 ,color=RED)
        cos_func.add_updater(lambda mob,dt: mob.become(FunctionGraph(lambda t: np.cos(t-dt) + 0.5, x_range =[-2*PI, 6*PI+self.renderer.time+dt] ,color=RED)))
        cos_func.add_updater(lambda mob,dt: mob.become(mob.shift(LEFT*(self.renderer.time+dt-1))))

        shifted_cos_func = cos_func.copy()
        shifted_cos_func.add_updater(lambda mob,dt: mob.become(FunctionGraph(lambda t: np.cos(t+(3*PI/2)-dt) + 0.5, x_range =[-2*PI, 6*PI+self.renderer.time+dt] ,color=RED)))
        shifted_cos_func.add_updater(lambda mob,dt: mob.become(mob.shift(LEFT*(self.renderer.time+dt-1))))

        shiftedback_cos_func = cos_func.copy()
        shiftedback_cos_func.add_updater(lambda mob,dt: mob.become(FunctionGraph(lambda t: np.cos(t+(6*PI/2)-dt) + 0.5, x_range =[-2*PI, 6*PI+self.renderer.time+dt] ,color=RED)))
        shiftedback_cos_func.add_updater(lambda mob,dt: mob.become(mob.shift(LEFT*(self.renderer.time+dt-1))))

        V = MathTex(r"V_m\cos(wt+\phi)").move_to(2*UP).scale(0.7)  ; V[0][0:11].set_color(RED)
        I = MathTex(r"I_m\cos(wt+\phi)").move_to(RIGHT*3+UP*2).scale(0.7)  ; I[0][0:11].set_color(BLUE)

        Text_Inphase = Tex(r'Voltage and current are in phase, as illustrated in the phasor diagram.').move_to(RIGHT*3+UP).scale(0.5) 
        Text_Outphase_VI = MathTex(r"\text{The voltage and current are } 90^\circ \text{ out of phase}.").move_to(RIGHT*3+UP).scale(0.5) 
        Text_Outphase_V = MathTex(r"\text{Specifically, the current lags the voltage by } 90^\circ.", color=RED).move_to(RIGHT*3 + UP*0.5).scale(0.5) 
        Text_Outphase_I = MathTex(r"\text{Specifically, the current leads the voltage by } 90^\circ.", color=BLUE).move_to(RIGHT*3 + UP*0.5).scale(0.5) 

        frameboxP = SurroundingRectangle(Text_Inphase, buff = .1)
        frameboxV = SurroundingRectangle(Text_Outphase_V, buff = .1)
        frameboxI = SurroundingRectangle(Text_Outphase_I, buff = .1)
    #Circuits    
        Resistor = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to[R=$R$] (0,2)
                (-2,2) to[short,i=$i$] (0,2)
                to[short,-o] (-2,2)
                (0,0) to[short, -o] (-2,0)
                to[short] (0,0)
		        (-2,2) to[open, v=$v$] (-2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0,); Resistor[0][10:13].set_color(RED); Resistor[0][4:6].set_color(BLUE)
        Resistor_Fasor = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to[R=$R$] (0,2)
                (-2,2) to[short,i=$I$] (0,2)
                to[short,-o] (-2,2)
                (0,0) to[short, -o] (-2,0)
                to[short] (0,0)
		        (-2,2) to[open, v=$V$] (-2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0).shift(LEFT*4+UP*1.4).scale(0.4); Resistor_Fasor[0][10:13].set_color(RED); Resistor_Fasor[0][4:6].set_color(BLUE)    
        Inductor = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to[L=$L$] (0,2)
                (-2,2) to[short,i=$i$] (0,2)
                to[short,-o] (-2,2)
                (0,0) to[short, -o] (-2,0)
                to[short] (0,0)
		        (-2,2) to[open, v=$v$] (-2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0).move_to(UP*0.2); Inductor[0][10:13].set_color(RED); Inductor[0][4:6].set_color(BLUE)
        Inductor_Fasor = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to[L=$L$] (0,2)
                (-2,2) to[short,i=$I$] (0,2)
                to[short,-o] (-2,2)
                (0,0) to[short, -o] (-2,0)
                to[short] (0,0)
		        (-2,2) to[open, v=$V$] (-2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0).shift(LEFT*4.5+UP*1.4).scale(0.4); Inductor_Fasor[0][10:13].set_color(RED); Inductor_Fasor[0][4:6].set_color(BLUE) 
        Capacitor = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to[C=$C$] (0,2)
                (-2,2) to[short,i=$i$] (0,2)
                to[short,-o] (-2,2)
                (0,0) to[short, -o] (-2,0)
                to[short] (0,0)
		        (-2,2) to[open, v=$v$] (-2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0).move_to(UP*0.2); Capacitor[0][10:13].set_color(RED); Capacitor[0][4:6].set_color(BLUE)
        Capacitor_Fasor = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to[C=$C$] (0,2)
                (-2,2) to[short,i=$I$] (0,2)
                to[short,-o] (-2,2)
                (0,0) to[short, -o] (-2,0)
                to[short] (0,0)
		        (-2,2) to[open, v=$V$] (-2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0).shift(LEFT*4.5+UP*1.4).scale(0.4); Capacitor_Fasor[0][10:13].set_color(RED); Capacitor_Fasor[0][4:6].set_color(BLUE)
        Resistor_VSource = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to[R=$R$] (0,2)
                (-2,2) to[short,i=$i$] (0,2)
                to[short,-o] (-2,2)
                (0,0) to[short, -o] (-2,0)
                to[short] (0,0)
		        (-2,2) to [sV] (-2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0,).shift(LEFT*4.15+UP*1.4).scale(0.4); Resistor_VSource[0][6:15].set_color(RED); Resistor_VSource[0][4:6].set_color(BLUE)
        Inductor_VSource = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to[L=$L$] (0,2)
                (-2,2) to[short,i=$i$] (0,2)
                to[short,-o] (-2,2)
                (0,0) to[short, -o] (-2,0)
                to[short] (0,0)
		        (-2,2) to [sV]  (-2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0).shift(LEFT*4.5+UP*1.4).scale(0.4); Inductor_VSource[0][10:13].set_color(RED); Inductor_VSource[0][4:6].set_color(BLUE)
        Capacitor_VSource = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to[C=$C$] (0,2)
                (-2,2) to[short,i=$i$] (0,2)
                to[short,-o] (-2,2)
                (0,0) to[short, -o] (-2,0)
                to[short] (0,0)
		        (-2,2) to [sV] (-2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0).shift(LEFT*4.5+UP*1.4).scale(0.4); Capacitor_VSource[0][10:13].set_color(RED); Capacitor_VSource[0][4:6].set_color(BLUE)
    #Equations below circuits
        Resistor_v = MathTex(
            r"v=iR", substrings_to_isolate="v,i"
            ).shift(DOWN*2.2+LEFT*0.5).scale(1.3).set_color_by_tex("v", RED).set_color_by_tex("i", BLUE) 
        Resistor_v_Fasor = MathTex(
            r"\overrightarrow{V}=\overrightarrow{I}R"
         ).shift(LEFT*4.5+UP*0.5).scale(0.5); Resistor_v_Fasor[0][0:3].set_color(RED); Resistor_v_Fasor[0][4:7].set_color(BLUE)
        Inductor_v = MathTex(
            r"v=L\frac{di}{dt}" #, substrings_to_isolate="x"
            ).shift(DOWN*2.4+LEFT*0.65).scale(1.3) ; Inductor_v[0][0:1].set_color(RED); Inductor_v[0][4:5].set_color(BLUE)
        Inductor_v_Fasor = MathTex(
            r"\overrightarrow{V}=jwL\overrightarrow{I}" #, substrings_to_isolate="x"
            ).shift(LEFT*4.5+UP*0.5).scale(0.5); Inductor_v_Fasor[0][0:3].set_color(RED); Inductor_v_Fasor[0][7:10].set_color(BLUE)
        Capacitor_v = MathTex(
            r"i=C\frac{dv}{dt}" #, substrings_to_isolate="x"
            ).shift(DOWN*2.4+LEFT*0.5).scale(1.3); Capacitor_v[0][0:1].set_color(BLUE); Capacitor_v[0][4:5].set_color(RED)
        Capacitor_v_Fasor = MathTex(
                r"\overrightarrow{I}=jwC\overrightarrow{V}" #, substrings_to_isolate="x"
            ).shift(LEFT*4.5+UP*0.5).scale(0.5); Capacitor_v_Fasor[0][0:3].set_color(BLUE); Capacitor_v_Fasor[0][7:10].set_color(RED)
    #Equations Resistor
        R_V_tempo = Resistor_v.copy().shift(LEFT*3.85+UP*2.7).scale(0.4)
        R_V_time_eq = MathTex(r"V_m\cos(wt+\phi)=I_m\cos(wt+\phi)*R").move_to(RIGHT*2+UP).scale(0.7)  ; R_V_time_eq[0][0:11].set_color(RED); R_V_time_eq[0][12:23].set_color(BLUE)
        R_V_phasor_eq = MathTex(r"V_m\angle\phi=I_m\angle\phi *R").move_to(RIGHT*2+UP*0.5).scale(0.7)
        R_V_vector_eq = MathTex(r"\overrightarrow{V}=\overrightarrow{I}R").move_to(RIGHT*2).scale(0.7); R_V_vector_eq[0][0:3].set_color(RED); R_V_vector_eq[0][4:7].set_color(BLUE)
    #Equations Inductor
        L_V_tempo = Inductor_v.copy().move_to(LEFT*4.5+UP*0.3).scale(0.4)
        L_V_time_eq = MathTex(r"V_m\cos(wt+\phi)=L\frac{d(I_m\cos(wt+\phi))}{dt}").move_to(RIGHT*2+UP*0.4).scale(0.7)  ; L_V_time_eq[0][0:11].set_color(RED); L_V_time_eq[0][15:26].set_color(BLUE)
        L_V_deriv_time_eq = MathTex(r"V_m\cos(wt+\phi)=-wL*I_m\sin(wt+\phi)").move_to(RIGHT*2+DOWN*0.3).scale(0.7)  ; L_V_deriv_time_eq[0][0:11].set_color(RED); L_V_deriv_time_eq[0][16:28].set_color(BLUE)
        L_V_deriv_trans_time_eq = MathTex(r"V_m\cos(wt+\phi)=wL*I_m\cos(wt+90+\phi)").move_to(RIGHT*2+DOWN*.3).scale(0.7)  ; L_V_deriv_trans_time_eq[0][0:11].set_color(RED); L_V_deriv_trans_time_eq[0][15:29].set_color(BLUE)
        L_V_phasor_eq = MathTex(r"V_m\angle\phi=wL*I_me^{j(90+\phi)}").move_to(RIGHT*2+DOWN*.3).scale(0.7)
        L_V_vector_eq = MathTex(r"\overrightarrow{V}=wL*I_me^{j90}e^{j\phi}").move_to(RIGHT*2+DOWN*.8).scale(0.7); L_V_vector_eq[0][0:3].set_color(RED); L_V_vector_eq[0][7:17].set_color(BLUE)
        #L_V_vector_p2_eq = MathTex(r"\overrightarrow{V}=wL*I_m\angle\phi e^{j90}").move_to(RIGHT*2+DOWN*.3).scale(0.7); L_V_vector_p2_eq[0][0:3].set_color(RED); L_V_vector_p2_eq[0][7:17].set_color(BLUE)
        L_V_vector_p3_eq = MathTex(r"\overrightarrow{V}=jwL\overrightarrow{I}").move_to(RIGHT*2+DOWN*.3).scale(0.7); L_V_vector_p3_eq[0][0:3].set_color(RED); L_V_vector_p3_eq[0][7:17].set_color(BLUE) ; L_V_vector_p3_eq[0][4:5].set_color(BLUE)              
    #Equations Capacitor
        C_V_tempo = Capacitor_v.copy().move_to(LEFT*4.5+UP*0.3).scale(0.4)
        C_V_time_eq = MathTex(r"I_m\cos(wt+\phi)=C\frac{d(V_m\cos(wt+\phi))}{dt}").move_to(RIGHT*2+UP*0.4).scale(0.7)  ; C_V_time_eq[0][0:11].set_color(BLUE); C_V_time_eq[0][15:26].set_color(RED)
        C_V_deriv_time_eq = MathTex(r"I_m\cos(wt+\phi)=-wC*V_m\sin(wt+\phi)").move_to(RIGHT*2+DOWN*0.3).scale(0.7)  ; C_V_deriv_time_eq[0][0:11].set_color(BLUE); C_V_deriv_time_eq[0][16:28].set_color(RED)
        C_V_deriv_trans_time_eq = MathTex(r"I_m\cos(wt+\phi)=wC*V_m\cos(wt+90+\phi)").move_to(RIGHT*2+DOWN*.3).scale(0.7)  ; C_V_deriv_trans_time_eq[0][0:11].set_color(BLUE); C_V_deriv_trans_time_eq[0][15:29].set_color(RED)
        C_V_phasor_eq = MathTex(r"I_m\angle\phi=wC*V_me^{j(90+\phi)}").move_to(RIGHT*2+DOWN*.3).scale(0.7)
        C_V_vector_eq = MathTex(r"\overrightarrow{I}=wC*V_me^{j90}e^{j\phi}").move_to(RIGHT*2+DOWN*.8).scale(0.7); C_V_vector_eq[0][0:3].set_color(BLUE); C_V_vector_eq[0][7:17].set_color(RED)
        #C_V_vector_p2_eq = MathTex(r"\overrightarrow{I}=wC*V_m\angle\phi e^{j90}").move_to(RIGHT*2+DOWN*.3).scale(0.7)
        C_V_vector_p3_eq = MathTex(r"\overrightarrow{I}=jwC\overrightarrow{V}").move_to(RIGHT*2+DOWN*.3).scale(0.7); C_V_vector_p3_eq[0][0:3].set_color(BLUE); C_V_vector_p3_eq[0][7:17].set_color(RED); C_V_vector_p3_eq[0][4:5].set_color(RED)   
    #Graph
        x_axis_R = Line(ORIGIN, RIGHT*3).add_tip()
        y_axis_R = Line(ORIGIN, UP*3).add_tip()
        vector_I_R = (Line(ORIGIN, 2*RIGHT).set_color(BLUE).add_tip()
                    .rotate(45*DEGREES, about_point=ORIGIN))
        vector_V_R = (Line(ORIGIN, 3*RIGHT).set_color(RED).add_tip()
                    .rotate(45*DEGREES, about_point=ORIGIN))          
        arc_phi_R = (Arc(radius=1, start_angle=0, angle=45*DEGREES)
                   .add_tip(tip_length=0.2, tip_width=0.2))
        #right_angle = RightAngle(vector_I, vector_V, length=0.3)
        labels_R = [
            #Tex("Re").next_to(RIGHT*3, DOWN),
            #Tex("Im").next_to(UP*3, LEFT),
            Tex("I").next_to(vector_I_R.get_end(), RIGHT),
            Tex("V").next_to(vector_V_R.get_end(), LEFT),
            Tex(r"$\phi$").next_to(arc_phi_R.point_from_proportion(.9), RIGHT)
            #Tex(r"$\omega$").next_to(arc_omega.point_from_proportion(.8), RIGHT)    
        ]
        graph_R = VGroup(x_axis_R, y_axis_R, *labels_R, vector_I_R, vector_V_R, arc_phi_R).shift(LEFT*4.3+DOWN*2.5).scale(0.4)
    #Circuits
        R_Symbol = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to [R=$R$] (2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0).scale(0.6).move_to(LEFT*3+DOWN); 
        L_Symbol = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to [L=$L$] (2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0).scale(0.6).move_to(DOWN); 
        C_Symbol = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to [C=$C$] (2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0).scale(0.6).move_to(RIGHT*3+DOWN); 
        Resistor = Tex(r"""
            \begin{circuitikz}[american]
            \ctikzset{capacitors/scale=0.6}
            \draw
		        (0,0) to[R=$R$] (0,2)
                (-2,2) to[short,i=$i$] (0,2)
                to[short,-o] (-2,2)
                (0,0) to[short, -o] (-2,0)
                to[short] (0,0)
		        (-2,2) to[open, v=$v$] (-2,0);
            \end{circuitikz}
            """
            ,stroke_width=2,tex_template = template,fill_opacity=0,); Resistor[0][10:13].set_color(RED); Resistor[0][4:6].set_color(BLUE)        
    #Resistor Animation    
        with self.voiceover(text="Welcome. In this video, we will learn an intriguing exploration, about the phasor relationships for circuit elements.") as tracker:
            self.play(Write(Title))
            self.play(Title.animate.shift(UP*3).scale(0.8))
        self.wait(2)
        with self.voiceover(text="Today, we'll uncover the fascinating interplay of voltage and current in the realm of electrical circuits.") as tracker:
            self.play(Create(VGroup(sin_func,cos_func)), run_time=3)
            self.play(Write(VGroup(V,I)))
            self.wait(2)
        with self.voiceover(text="Much like the symbiotic relationships found in nature.") as tracker:
            self.play(FadeOut(VGroup(V,I)))
        with self.voiceover(text="Circuit elements such as resistors.") as tracker:
            self.play(Write(R_Symbol))
            self.play(ReplacementTransform(cos_func,shifted_cos_func))
        with self.voiceover(text="Inductors.") as tracker:
            self.play(Write(L_Symbol))
            self.play(ReplacementTransform(shifted_cos_func,cos_func))
        with self.voiceover(text="And Capacitors, interact with alternating currents in a mesmerizing choreography of phases.") as tracker:
            self.play(Write(C_Symbol))
            self.play(ReplacementTransform(cos_func,shiftedback_cos_func))
            self.play(FadeOut(R_Symbol,L_Symbol,C_Symbol))
        self.wait(3)
        with self.voiceover(text="Our journey begins with the humble resistor, where voltage and current maintain a synchronous relationship, moving in lockstep.") as tracker:
            self.play(Write(VGroup(Resistor,Resistor_v)))
            self.play(ReplacementTransform(shiftedback_cos_func,shifted_cos_func))
            self.wait(2)
            self.play(Write(Resistor_Title))
            self.play(Create(frameboxRES))
            self.wait(2)
            self.play(FadeOut(Resistor_Title,frameboxRES))
        self.wait(1.5)
        with self.voiceover(text="Imagine a resistor, that unassuming yet vital component of an electrical circuit, connected to an AC voltage source.") as tracker:
            sin_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1/4)))
            shifted_cos_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1/4)))
            self.play(VGroup(Resistor,Resistor_v).animate.shift(LEFT*4+UP*1.6).scale(0.4))
            self.play(R_V_tempo.animate.move_to(RIGHT*2+UP*1.5).scale(1.3))
            self.play(Transform(Resistor,Resistor_VSource))
        self.wait(1.5)
        with self.voiceover(text="Ohm's Law, that essential principle of electrical circuits, reigns supreme here, dictating that voltage is directly proportional to current for a given resistance.") as tracker:
            self.camera.frame.save_state()
            self.play(self.camera.frame.animate.set(width=R_V_tempo.width*6).move_to(R_V_tempo).shift(DOWN))
        self.wait(1.5)   
        with self.voiceover(text="In the time domain, we deal with the instantaneous values of voltage and current, which vary continuously as a function of time.") as tracker:
            self.play(Write(R_V_time_eq))
        self.wait(1.5)   
        with self.voiceover(text="Here, the phasor domain comes to our aid, transforming those time-dependent sine or cosine functions into a more manageable form: the phasor.") as tracker:
            self.play(Write(R_V_phasor_eq))
            self.wait(3)
            self.play(Write(R_V_vector_eq))
            self.play(FadeOut(VGroup(R_V_phasor_eq,R_V_time_eq,R_V_tempo)))
            self.play(R_V_vector_eq.animate.shift(UP*1.5))
            self.wait(3)
            self.play(Restore(self.camera.frame))
        with self.voiceover(text="For resistors, the beauty lies in the simplicity of their behavior. Both the voltage and current maintain the same phase in the phasor domain.") as tracker:
            self.add(VGroup(Resistor_Fasor,Resistor_v_Fasor))
            self.play(Resistor_v_Fasor.animate.shift(RIGHT*2), Resistor_Fasor.animate.shift(RIGHT*1.5))
            self.play(R_V_vector_eq.animate.shift(LEFT*4.5+DOWN).scale(0.7))  
        with self.voiceover(text="This phase alignment means that the phasor representation of voltage, V, and current, I, are collinear in the phasor diagram, preserving their direct relationship.") as tracker:
            self.play(Write(graph_R))
            self.bring_to_front(vector_I_R)
            self.play(Write(Text_Inphase))
            self.play(Create(frameboxP))
            self.wait(3)
        with self.voiceover(text="Thus, Ohm's Law maintains its integrity even in the phasor domain.") as tracker:
            self.play(FadeOut(VGroup(Resistor,Resistor_v,Resistor_Fasor,Resistor_v_Fasor,R_V_vector_eq,graph_R,Text_Inphase,frameboxP)))
            sin_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1)))
            cos_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1)))
        self.wait(3)
    #Inductor Animation
        with self.voiceover(text="Ah, the inductor - that curious component that introduces the element of time into electrical circuits. Let's explore, the intriguing relationship between, the time domain and phasor domain for inductors.") as tracker:
            self.play(Write(VGroup(Inductor,Inductor_v)))
            self.play(ReplacementTransform(shifted_cos_func,cos_func))
            self.play(Write(Inductor_Title))
            self.play(Create(frameboxIND))
            self.wait(3)
            self.play(FadeOut(Inductor_Title,frameboxIND))
        self.wait(1.5)
        with self.voiceover(text="In the time domain, the behavior of an inductor is governed by Faraday's Law of Induction. This law dictates that, the voltage across an inductor is proportional to the rate of change of current through it.") as tracker:
            sin_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1/4)))
            cos_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1/4)))
            self.play(VGroup(Inductor,Inductor_v).animate.shift(LEFT*4.2+UP*1.6).scale(0.4))
            self.play(L_V_tempo.animate.move_to(RIGHT*2+UP*1.5).scale(1.3))    
        self.wait(1.5)
        with self.voiceover(text="However, when dealing with AC circuits, the time domain expressions become more complex due to the continuously changing nature of alternating currents.") as tracker:
            self.play(Transform(Inductor,Inductor_VSource))
            self.play(Write(L_V_time_eq))
            self.camera.frame.save_state()
            self.play(self.camera.frame.animate.set(width=L_V_tempo.width*6).move_to(L_V_tempo).shift(DOWN))
            self.wait(0.3)
            self.play(Write(L_V_deriv_time_eq))
            self.play(FadeOut(L_V_time_eq))
        self.wait(1.5)
        with self.voiceover(text="After some algeabric manipulation. We can calculate the differentiation of the current value to a more standard form.") as tracker:
            self.play(L_V_deriv_time_eq.animate.shift(UP*0.5))
            self.wait(1)
            self.play(Write(L_V_deriv_trans_time_eq))
            self.play(FadeOut(L_V_deriv_time_eq))
            self.play(L_V_deriv_trans_time_eq.animate.shift(UP*0.5))
        self.wait(1.5)
        with self.voiceover(text="This is where the phasor domain steps in, simplifying the analysis.") as tracker:
            self.play(Write(L_V_phasor_eq))
            self.wait(1)
            self.play(ReplacementTransform(L_V_deriv_time_eq,L_V_vector_eq))
            self.play(FadeOut(L_V_phasor_eq))
            self.play(L_V_vector_eq.animate.shift(UP*0.5))
        self.wait(1.5)    
        with self.voiceover(text="To transition between the time and frequency domains, Euler's formula plays a crucial role in the realm of signal processing and Fourier analysis.") as tracker:
            self.play(FadeOut(L_V_deriv_trans_time_eq))
            self.play(L_V_vector_eq.animate.shift(UP*0.5))
            self.wait(1)
            self.play(Write(L_V_vector_p3_eq))
        self.wait(1.5)
        with self.voiceover(text="For an inductor, the current lags behind the voltage by 90 degrees in the phasor domain, creating a fascinating contrast to the in-phase relationship of resistors. In other words, when the voltage is at its peak, the current is just beginning its journey.") as tracker:
            self.play(FadeOut(L_V_tempo,L_V_vector_eq))
            self.play(L_V_vector_p3_eq.animate.shift(UP*1.8))
            self.wait(1)
            self.play(Restore(self.camera.frame))
            self.add(VGroup(Inductor_Fasor,Inductor_v_Fasor))
            self.play(VGroup(Inductor_Fasor,Inductor_v_Fasor).animate.shift(RIGHT*1.5))
            self.play(L_V_vector_p3_eq.animate.shift(LEFT*5+DOWN).scale(0.7))
        self.wait(1.5)
        with self.voiceover(text="This phase difference is captured elegantly in the phasor diagram, where the current phasor lags the voltage phasor by a quarter of a cycle. Thus, the phasor domain beautifully encapsulates the time-dependent behavior of inductors.") as tracker:
            self.play(Write(graph_R))
            self.play(labels_R[1].animate.shift(LEFT*1.7),Rotate(vector_V_R, PI/2, about_point=x_axis_R.get_start()))
            self.play(Write(Text_Outphase_VI))
            self.play(Write(Text_Outphase_V))
            self.play(Create(frameboxV))
            self.wait(3)
            self.play(FadeOut(VGroup(Inductor,Inductor_v,Inductor_Fasor,Inductor_v_Fasor,L_V_vector_p3_eq,graph_R,Text_Outphase_VI,Text_Outphase_V,frameboxV)))
            self.wait(3) 
            sin_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1)))
            cos_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1)))
            shifted_cos_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1)))
            self.play(ReplacementTransform(cos_func, shifted_cos_func))
        self.wait(3)    
    #Capacitor Animation 
        with self.voiceover(text="Let's now shed light on the intricate relationship between the time domain and phasor domain for capacitors, those essential components that store electrical energy.") as tracker:
            self.play(Write(VGroup(Capacitor,Capacitor_v)))
            self.play(ReplacementTransform(shifted_cos_func, shiftedback_cos_func))
            self.play(Write(Capacitor_Title))
            self.play(Create(frameboxCAP))
            self.wait(3)
            self.play(FadeOut(Capacitor_Title,frameboxCAP))
            sin_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1/4)))
            cos_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1/4)))
            shifted_cos_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1/4)))
            shiftedback_cos_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1/4)))
        self.wait(1.5)
        with self.voiceover(text="In the time domain, a capacitor's behavior is dictated by the relationship between voltage and the rate of change of charge. The capacitor current, is proportional to the rate of change of voltage across its terminals.") as tracker:
            self.play(VGroup(Capacitor,Capacitor_v).animate.shift(LEFT*4.2+UP*1.6).scale(0.4))
            self.play(C_V_tempo.animate.move_to(RIGHT*2+UP*1.5).scale(1.3))
            self.play(Transform(Capacitor,Capacitor_VSource))
            self.play(Write(C_V_time_eq))
        self.wait(1.5)    
        with self.voiceover(text="However, when we enter the realm of AC circuits, the time domain expressions become more complex, necessitating the use of the phasor domain.") as tracker:
            self.camera.frame.save_state()
            self.play(self.camera.frame.animate.set(width=C_V_tempo.width*6).move_to(C_V_tempo).shift(DOWN))
            self.wait(0.3)
            self.play(Write(C_V_deriv_time_eq))
            self.play(FadeOut(C_V_time_eq))
        self.wait(1.5)  
        with self.voiceover(text="Similary to the inductor case, after some algeabric manipulation. We can calculate the differentiation of the voltage value to a more standard form.") as tracker:
            self.play(C_V_deriv_time_eq.animate.shift(UP*0.5))
            self.wait(1)
            self.play(Write(C_V_deriv_trans_time_eq))
            self.play(FadeOut(C_V_deriv_time_eq))
            self.play(C_V_deriv_trans_time_eq.animate.shift(UP*0.5))
        self.wait(1.5)
        with self.voiceover(text="This is where the phasor domain steps in, simplifying the analysis.") as tracker:
            self.play(Write(C_V_phasor_eq))
            self.wait(1)
            self.play(ReplacementTransform(C_V_deriv_time_eq,C_V_vector_eq))
        self.wait(1.5)
        with self.voiceover(text="As mentioned before, Euler's formula plays a crucial role in this simplification.") as tracker:
            self.play(FadeOut(C_V_phasor_eq))
            self.play(C_V_vector_eq.animate.shift(UP*0.5))
            self.play(FadeOut(C_V_deriv_trans_time_eq))
            self.play(C_V_vector_eq.animate.shift(UP*0.5))
            self.wait(1)
            self.play(Write(C_V_vector_p3_eq))
        self.wait(1.5)
        with self.voiceover(text="For a capacitor, the current leads the voltage by 90 degrees in the phasor domain. This intriguing phase relationship implies that when the current through the capacitor is at its peak, the voltage is just beginning to build up.") as tracker:
            self.wait(1)
            self.play(FadeOut(C_V_vector_eq,C_V_tempo))
            self.play(C_V_vector_p3_eq.animate.shift(UP*1.8))
            self.wait(1)
            self.play(Restore(self.camera.frame))
            self.add(VGroup(Capacitor_Fasor,Capacitor_v_Fasor))
            self.play(VGroup(Capacitor_Fasor,Capacitor_v_Fasor).animate.shift(RIGHT*1.5))
            self.play(C_V_vector_p3_eq.animate.shift(LEFT*5+DOWN).scale(0.7))
        self.wait(1.5)
        with self.voiceover(text="The phasor diagram captures this lead-lag relationship visually, with the current phasor leading the voltage phasor by a quarter of a cycle. In summary, the phasor domain encapsulates the time-dependent behavior of capacitors, highlighting the essential unity of electrical principles across domains.") as tracker:
            self.play(Write(graph_R))
            self.play(labels_R[1].animate.shift(RIGHT*1.7), Rotate(vector_V_R, -PI/2, about_point=x_axis_R.get_start()))
            self.play(labels_R[0].animate.shift(LEFT*1.4), Rotate(vector_I_R, PI/2, about_point=x_axis_R.get_start()))
            self.play(Write(Text_Outphase_VI))
            self.play(Write(Text_Outphase_I))
            self.play(Create(frameboxI))
            self.wait(3)
        self.wait(1.5)
        with self.voiceover(text="That's all for now. Thank you for watching. Don't forget to leave a like and subscribe for more videos like this.") as tracker:
            sin_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1)))
            cos_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1)))
            shifted_cos_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1)))
            shiftedback_cos_func.add_updater(lambda mob,dt: mob.become(mob.set_stroke(opacity=1)))
            self.play(FadeOut(VGroup(Capacitor,Capacitor_v,Capacitor_Fasor,Capacitor_v_Fasor,C_V_vector_p3_eq,graph_R,Text_Outphase_VI,Text_Outphase_I,frameboxI)))
        self.wait(1.5)

            