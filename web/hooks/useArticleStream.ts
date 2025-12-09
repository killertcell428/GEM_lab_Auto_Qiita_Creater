'use client';

import { useEffect, useState, useRef } from 'react';

interface StreamEvent {
  type: string;
  phase?: string;
  status_text?: string;
  message?: string;
  article_id?: string;
}

export function useArticleStream(articleId: string | null) {
  const [events, setEvents] = useState<StreamEvent[]>([]);
  const [connected, setConnected] = useState(false);
  const eventSourceRef = useRef<EventSource | null>(null);

  useEffect(() => {
    if (!articleId) {
      return;
    }

    const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';
    const eventSource = new EventSource(`${apiUrl}/api/articles/${articleId}/stream`);
    eventSourceRef.current = eventSource;

    eventSource.onopen = () => {
      setConnected(true);
    };

    eventSource.onmessage = (event) => {
      try {
        // SSEイベントのdata:プレフィックスを除去
        let dataStr = event.data;
        if (dataStr.startsWith('data: ')) {
          dataStr = dataStr.substring(6); // 'data: 'を除去
        }
        const data: StreamEvent = JSON.parse(dataStr);
        setEvents((prev) => [...prev, data]);
      } catch (err) {
        console.error('Failed to parse SSE event:', err, 'Raw data:', event.data);
      }
    };

    eventSource.onerror = () => {
      setConnected(false);
      eventSource.close();
    };

    return () => {
      eventSource.close();
    };
  }, [articleId]);

  return { events, connected };
}

