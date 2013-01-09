classdef Animation < handle
  properties
    figure;
    paused;
    timer;
    frame;
    length;
    render;
    fps;
  end

  methods
    function this = Animation(fig)
      this.figure = fig;
      this.paused = true;
      this.timer = [];
      this.frame = 1;
      this.length = 0;
      this.fps = 25;
      this.render = [];

      set(this.figure, 'KeyPressFcn', @(obj, evt) this.key_press(obj, evt), ...
        'CloseRequestFcn', @(src, evt) this.request_close_figure(src, evt));
    end

    function key_press(this, obj, evt)
      if strcmp(evt.Key, 'space')
        if this.paused
          this.play();
        else
          this.pause();
        end
      elseif strcmp(evt.Key, 'leftarrow')
        this.pause();
        this.prev_frame();
        this.render_frame();
      elseif strcmp(evt.Key, 'rightarrow')
        this.pause();
        this.next_frame();
        this.render_frame();
      end
    end

    function request_close_figure(this, src, evt)
      if isvalid(this.timer)
        stop(this.timer);
        delete(this.timer);
      end
      delete(this.figure);
    end

    function play(this)
      this.paused = false;
      this.timer = timer('TimerFcn', @(x, y) this.timer_fire(), ...
        'Period', round(1000 / this.fps) / 1000, 'ExecutionMode', 'fixedRate');

      this.timer_fire();
      start(this.timer);
    end

    function pause(this)
      this.paused = true;

      % Stop and delete timer if it exists.
      if isvalid(this.timer)
        stop(this.timer);
        delete(this.timer);
      end
    end

    function timer_fire(this)
      this.next_frame();
      this.render_frame();
    end

    function render_frame(this)
      this.render(this.figure, this.frame);
      drawnow;
    end

    function advance_frame(this, dt)
      this.frame = mod(this.frame + dt - 1, this.length) + 1;
    end

    function next_frame(this)
      this.advance_frame(1);
    end

    function prev_frame(this)
      this.advance_frame(-1);
    end
  end
end
