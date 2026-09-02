import { useMemo, useState } from "react";
import { Link } from "react-router-dom";

import Icon from "../components/Icon";
import {
  KNOWLEDGE_CATEGORIES,
  searchKnowledge,
  type KnowledgeCategory,
} from "../data/knowledge";

type CategoryFilter = "All" | KnowledgeCategory;

export default function Learn() {
  const [query, setQuery] = useState("");
  const [category, setCategory] = useState<CategoryFilter>("All");
  const topics = useMemo(() => searchKnowledge(query, category), [category, query]);

  const clear = () => {
    setQuery("");
    setCategory("All");
  };

  return (
    <>
      <div className="topbar">
        <div className="grow">
          <h1>Learn</h1>
          <p className="sub">Curated molecular biology and genetic engineering knowledge, available offline.</p>
        </div>
      </div>

      <div className="content learn-library">
        <section className="card learn-intro" aria-labelledby="learn-title">
          <div className="card-body">
            <div className="learn-intro-copy">
              <span className="label">Deterministic reference</span>
              <h2 id="learn-title">Ask the library, not a black box</h2>
              <p>
                Search practical explanations used by G-Synth. The content is fixed, reviewable and
                restricted to molecular biology and genetic engineering; no sequence or question is
                sent to an external AI service.
              </p>
            </div>
            <div className="learn-scope" aria-label="Knowledge scope">
              <span><Icon name="check" size={15} /> 10 reviewed topics</span>
              <span><Icon name="check" size={15} /> Experimental limits stated</span>
              <span><Icon name="check" size={15} /> Module links included</span>
            </div>
          </div>
        </section>

        <section className="learn-search" aria-label="Search the knowledge library">
          <div className="field learn-query">
            <label htmlFor="learn-query">Search by question or keyword</label>
            <div className="learn-query-control">
              <Icon name="search" size={18} />
              <input
                id="learn-query"
                type="text"
                value={query}
                onChange={(event) => setQuery(event.target.value)}
                placeholder="e.g. Why does a 5′ primer tail not hybridise in cycle 1?"
              />
              {query && <button className="btn btn-ghost" onClick={() => setQuery("")}>Clear search</button>}
            </div>
          </div>
          <div className="learn-categories" role="group" aria-label="Filter by subject">
            {KNOWLEDGE_CATEGORIES.map((item) => (
              <button
                key={item}
                className={`btn ${category === item ? "btn-primary" : "btn-outline"}`}
                aria-pressed={category === item}
                onClick={() => setCategory(item)}
              >
                {item}
              </button>
            ))}
          </div>
        </section>

        <div className="learn-results-head" role="status" aria-live="polite">
          <strong>{topics.length} {topics.length === 1 ? "topic" : "topics"}</strong>
          {(query || category !== "All") && <button className="btn btn-ghost" onClick={clear}>Reset filters</button>}
        </div>

        {topics.length ? (
          <div className="learn-topic-grid">
            {topics.map((topic, index) => (
              <details className="card learn-topic" key={topic.id} open={Boolean(query) && index === 0}>
                <summary>
                  <span className="learn-topic-index">{String(index + 1).padStart(2, "0")}</span>
                  <span className="learn-topic-copy">
                    <span className="label">{topic.category}</span>
                    <strong>{topic.title}</strong>
                    <small>{topic.summary}</small>
                  </span>
                  <Icon name="chevronDown" size={18} />
                </summary>
                <div className="learn-topic-body">
                  <h3>Essential points</h3>
                  <ul>
                    {topic.essentials.map((point) => <li key={point}>{point}</li>)}
                  </ul>
                  <div className="notice notice-info learn-caution">
                    <strong>Interpretation limit.</strong> {topic.caution}
                  </div>
                  <div className="learn-topic-actions">
                    <Link className="btn btn-primary" to={topic.route}>{topic.routeLabel}</Link>
                    {topic.references.map((reference) => (
                      <a key={reference.url} href={reference.url} target="_blank" rel="noreferrer">
                        {reference.label}
                      </a>
                    ))}
                  </div>
                </div>
              </details>
            ))}
          </div>
        ) : (
          <div className="card empty learn-empty">
            <Icon name="search" size={34} />
            <h2>No reviewed topic matches</h2>
            <p>Try a shorter molecular-biology term such as “primer”, “HindIII”, “gel”, “annotation” or “sequencing”.</p>
            <button className="btn btn-outline" onClick={clear}>Show all topics</button>
          </div>
        )}
      </div>
    </>
  );
}
