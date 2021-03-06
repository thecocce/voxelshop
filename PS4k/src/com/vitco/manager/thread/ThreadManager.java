package com.vitco.manager.thread;

import com.vitco.manager.action.ActionManager;
import com.vitco.manager.action.ChangeListener;
import com.vitco.manager.action.types.StateActionPrototype;
import org.springframework.beans.factory.annotation.Autowired;

import javax.annotation.PostConstruct;
import java.util.ArrayList;

/**
 * Helps with creating threads.
 */
public class ThreadManager implements ThreadManagerInterface {

    private ActionManager actionManager;
    // set the action handler
    @Override
    @Autowired
    public final void setActionManager(ActionManager actionManager) {
        this.actionManager = actionManager;
    }

    private final ArrayList<LifeTimeThread> threads = new ArrayList<LifeTimeThread>();

    // manage a thread
    @Override
    public synchronized void manage(LifeTimeThread thread) {
        threads.add(thread);
        thread.start();
    }

    // remove a thread from the managed list
    // (this also stops the execution of this thread)
    @Override
    public synchronized void remove(LifeTimeThread thread) {
        thread.stopThread();
        threads.remove(thread);
    }

    private final ThreadManager thisInstance = this;

    @PostConstruct
    @Override
    public void init() {
        // close all threads when the program exits
        actionManager.performWhenActionIsReady("program_closing_event", new Runnable() {
            @Override
            public void run() {
                ((StateActionPrototype)actionManager.getAction("program_closing_event")).addChangeListener(new ChangeListener() {
                    @Override
                    public void actionFired(boolean b) {
                        if (b) {
                            synchronized (thisInstance) {
                                for (LifeTimeThread thread : threads) {
                                    thread.stopThread();
                                }
                            }
                        }
                    }
                });
            }
        });
    }
}
