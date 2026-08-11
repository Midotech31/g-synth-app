import { useCallback, useEffect, useRef, useState } from "react";

import { ApiError, api } from "../api/client";
import Icon from "../components/Icon";

type Message = { role: "user" | "assistant"; content: string };

export default function Learn() {
  const [messages, setMessages] = useState<Message[]>([]);
  const [input, setInput] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");
  const bottomRef = useRef<HTMLDivElement>(null);
  const inputRef = useRef<HTMLTextAreaElement>(null);

  useEffect(() => {
    bottomRef.current?.scrollIntoView({ behavior: "smooth" });
  }, [messages, loading]);

  const send = useCallback(async () => {
    const question = input.trim();
    if (!question || loading) return;

    setInput("");
    setError("");
    const userMsg: Message = { role: "user", content: question };
    setMessages((prev) => [...prev, userMsg]);
    setLoading(true);

    try {
      const history = messages.map((m) => ({ role: m.role, content: m.content }));
      const { answer } = await api.askTutor(question, history);
      setMessages((prev) => [...prev, { role: "assistant", content: answer }]);
    } catch (err) {
      const detail =
        err instanceof ApiError
          ? err.message
          : "Something went wrong. Please try again.";
      setError(detail);
    } finally {
      setLoading(false);
      inputRef.current?.focus();
    }
  }, [input, loading, messages]);

  const handleKey = (e: React.KeyboardEvent) => {
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      void send();
    }
  };

  const clear = () => {
    setMessages([]);
    setError("");
    inputRef.current?.focus();
  };

  return (
    <>
      <div className="topbar">
        <div className="grow">
          <h1>Learn</h1>
          <p className="sub">
            Ask about gene synthesis, cloning, molecular biology, or genetic
            engineering.
          </p>
        </div>
        {messages.length > 0 && (
          <button className="btn btn-outline" onClick={clear}>
            New conversation
          </button>
        )}
      </div>

      <div className="chat-area">
        {/* A transcript that grows at the bottom is what `log` describes, and
            it is announced without interrupting whatever is being read. */}
        <div
          className="chat-messages"
          role="log"
          aria-live="polite"
          aria-relevant="additions"
          aria-busy={loading}
          aria-label="Conversation"
        >
          {messages.length === 0 && !loading && (
            <div className="chat-empty">
              <Icon name="book" size={40} />
              <h2>Study assistant</h2>
              <p className="note">
                Ask a question about gene synthesis, cloning, codon optimisation,
                restriction enzymes, or anything else in molecular biology and
                genetic engineering. The assistant is scoped to these topics and
                will point you to the right part of G-Synth when relevant.
              </p>
              <div className="chat-suggestions">
                {[
                  "What is a sticky end and why does it matter for cloning?",
                  "Explain codon optimisation in simple terms.",
                  "How does Merzoug assembly differ from Gibson assembly?",
                ].map((q) => (
                  <button
                    key={q}
                    className="btn btn-outline chat-suggestion"
                    onClick={() => {
                      setInput(q);
                      inputRef.current?.focus();
                    }}
                  >
                    {q}
                  </button>
                ))}
              </div>
            </div>
          )}

          {messages.map((msg, i) => (
            <div key={i} className={`chat-bubble chat-${msg.role}`}>
              <div className="chat-role label">
                {msg.role === "user" ? "You" : "Assistant"}
              </div>
              <div className="chat-text">{msg.content}</div>
            </div>
          ))}

          {loading && (
            <div className="chat-bubble chat-assistant">
              <div className="chat-role label">Assistant</div>
              <div className="chat-text chat-thinking">
                <span className="spinner" /> Thinking&hellip;
              </div>
            </div>
          )}

          {error && (
            <div className="notice notice-error" style={{ margin: "0.5rem 0" }} role="alert">
              {error}
            </div>
          )}

          <div ref={bottomRef} />
        </div>

        <div className="chat-input-bar">
          <textarea
            ref={inputRef}
            className="chat-input"
            placeholder="Ask a question…"
            aria-label="Your question"
            value={input}
            onChange={(e) => setInput(e.target.value)}
            onKeyDown={handleKey}
            rows={1}
            maxLength={2000}
            disabled={loading}
          />
          <button
            className="btn btn-primary chat-send"
            onClick={() => void send()}
            disabled={!input.trim() || loading}
            aria-label="Send question"
          >
            <Icon name="arrowRight" size={18} />
          </button>
        </div>
      </div>
    </>
  );
}
